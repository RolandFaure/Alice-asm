#include "reduce_and_expand.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <chrono>
#include <algorithm>
#include <omp.h> //for efficient parallelization
#include <set>
#include <unordered_set>
#include <map>
#include <atomic>
#include <memory>
#include <shared_mutex>
#include <mutex>
#include <zlib.h>

#include "robin_hood.h"
#include "basic_graph_manipulation.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

using std::cout;
using std::endl;
using std::string;
using std::vector;
using robin_hood::unordered_map;
using std::pair;
using std::cerr;
using std::ifstream;
using std::ofstream;

static void print_timestamp(){
    time_t now = time(0);
    tm *ltm = localtime(&now);
    cout << "[" << ltm->tm_mday << "/" << 1 + ltm->tm_mon << "/" << 1900 + ltm->tm_year << " " << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec << "]";
}

static inline bool is_ACGT(char c){
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

/**
 * @brief The parsers of the pipeline read single-line, uppercase FASTA files, and split them in chunks by seeking in
 * the file. Check if the input already is in that format, and if not (gzipped, FASTQ, multi-line or lowercase), write
 * a normalized copy in the tmp folder.
 *
 * @return path to the file to use for the rest of the pipeline
 */
string prepare_reads(string input_file, string tmp_folder){

    //is the file gzipped or a FASTQ?
    bool needs_normalization = false;
    {
        ifstream input(input_file, std::ios::binary);
        if (!input.is_open()){
            cerr << "ERROR: could not open file " << input_file << "\n";
            exit(1);
        }
        unsigned char magic[2] = {0, 0};
        input.read((char*) magic, 2);
        if (input.gcount() == 2 && magic[0] == 0x1f && magic[1] == 0x8b){
            needs_normalization = true;
        }
        else if (input.gcount() >= 1 && magic[0] != '>'){
            needs_normalization = true;
        }
    }

    //plain FASTA: check that each record is one header line followed by exactly one uppercase sequence line
    if (!needs_normalization){
        ifstream input(input_file);
        string line;
        bool previous_was_header = false;
        while (std::getline(input, line)){
            if (!line.empty() && line[0] == '>'){
                previous_was_header = true;
                continue;
            }
            if (!previous_was_header){ //sequence spread over several lines
                needs_normalization = true;
                break;
            }
            previous_was_header = false;
            for (char c : line){
                if (c < 'A' || c > 'Z'){ //lowercase, carriage return...
                    needs_normalization = true;
                    break;
                }
            }
            if (needs_normalization){
                break;
            }
        }
    }

    if (!needs_normalization){
        return input_file;
    }

    string normalized_file = tmp_folder + "reads.fa";
    cout << " - Converting the input reads to single-line uppercase FASTA in " << normalized_file << endl;
    gzFile fp = gzopen(input_file.c_str(), "r");
    if (fp == NULL){
        cerr << "ERROR: could not open file " << input_file << "\n";
        exit(1);
    }
    ofstream out(normalized_file);
    kseq_t *seq = kseq_init(fp);
    int l;
    while ((l = kseq_read(seq)) >= 0){
        if (seq->seq.l == 0){
            continue;
        }
        out << ">" << seq->name.s;
        if (seq->comment.l > 0){
            out << " " << seq->comment.s;
        }
        out << "\n";
        for (size_t i = 0 ; i < seq->seq.l ; i++){
            seq->seq.s[i] = toupper(seq->seq.s[i]);
        }
        out.write(seq->seq.s, seq->seq.l);
        out << "\n";
    }
    if (l < -1){
        cerr << "ERROR: the input file " << input_file << " is truncated or malformed\n";
        exit(1);
    }
    kseq_destroy(seq);
    gzclose(fp);
    out.close();
    if (!out){
        cerr << "ERROR: could not write " << normalized_file << "\n";
        exit(1);
    }
    return normalized_file;
}

/**
 * @brief Go through the reads of one chunk of a single-line FASTA file and call process(name_line, sequence) on each.
 * A read belongs to the chunk in which its header line starts.
 */
template <typename F>
static void for_each_read_in_chunk(const string& file, int chunk, unsigned long long size_of_chunk, F process){
    ifstream input(file);
    input.seekg(chunk*size_of_chunk);
    string line;
    string name_line;
    bool next_line_is_seq = false;
    while (std::getline(input, line)){
        if (line[0] == '>'){
            name_line = line;
            next_line_is_seq = true;
        }
        else if (next_line_is_seq){
            next_line_is_seq = false;
            process(name_line, line);
            long long position = input.tellg();
            if (position < 0 || (unsigned long long) position >= (chunk+1)*size_of_chunk){
                break;
            }
        }
    }
}

/**
 * @brief MSR-compress one read: roll a hash over the windows of `order` (homopolymer-compressed) bases and output one
 * base each time the canonical hash is a multiple of `compression`. Both the compression of the reads and the second
 * pass over the reads during the inflation must use this function, so that they sample exactly the same positions.
 *
 * @param on_sample called as on_sample(compressed_base, position of the middle of the window in the read)
 * @param on_break called when a non-ACGT base interrupts the read: the windows overlapping it are not sampled
 */
template <typename F, typename G>
static void sample_read(string& line, int order, int compression, bool homopolymer_compression, bool track_middle, F on_sample, G on_break){
    uint64_t hash_foward = 0;
    uint64_t hash_reverse = 0;
    size_t pos_end = 0;
    long pos_begin = -order;
    long pos_middle = track_middle ? -(order+1)/2 : -6666; //-6666 is the special value to not compute the middle position
    int number_of_valid_bases = 0; //number of consecutive ACGT (homopolymer-compressed) bases hashed

    while (pos_end < line.size()){
        char hashed_base = line[pos_end];
        roll(hash_foward, hash_reverse, order, line, pos_end, pos_begin, pos_middle, homopolymer_compression);
        if (!is_ACGT(hashed_base)){
            number_of_valid_bases = 0;
            on_break();
            continue;
        }
        number_of_valid_bases++;
        if (number_of_valid_bases > order){
            if (hash_foward < hash_reverse && hash_foward % compression == 0){
                on_sample("ACGT"[(hash_foward/compression)%4], pos_middle);
            }
            else if (hash_foward >= hash_reverse && hash_reverse % compression == 0){
                on_sample("TGCA"[(hash_reverse/compression)%4], pos_middle);
            }
        }
    }
}

/**
 * @brief MSR the input sequencing
 *
 * @param input_file single-line FASTA file (see prepare_reads)
 * @param output_file
 * @param order
 * @param compression
 * @param num_threads
 * @param homopolymer_compression
 **/
void reduce(string input_file, string output_file, int order, int compression, int num_threads, bool homopolymer_compression) {

    print_timestamp();
    cout << " Starting pipeline" << endl;

    std::ifstream input(input_file,std::ios::binary | std::ios::ate);
    if (!input.is_open())
    {
        std::cout << "Could not open file " << input_file << std::endl;
        exit(1);
    }

    std::streamoff file_size = input.tellg();
    unsigned long long size_of_chunk = 10000000;
    input.close();

    //clear the output file
    std::ofstream out(output_file);
    out.close();

    std::atomic<unsigned long long> seq_num (0);
    std::atomic<unsigned long long> output_limit (0);
    //parallelize on num_threads threads
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int chunk = 0 ; chunk <= file_size/size_of_chunk ; chunk++){

        string output_file_chunk = output_file + "_"+ std::to_string(chunk);
        std::ofstream out(output_file_chunk);

        for_each_read_in_chunk(input_file, chunk, size_of_chunk, [&](const string& name, string& line){
            unsigned long long n = seq_num++;
            if (n >= output_limit && omp_get_thread_num() == 0){
                #pragma omp critical
                {
                    print_timestamp();
                    cout << " Compressed " << n << " reads" << endl;
                }
                output_limit += 50000;
            }

            string name_line = name;
            bool first_base = true; //to output the name line when the first base is outputted (not before to avoid empty lines)
            sample_read(line, order, compression, homopolymer_compression, false,
                [&](char base, long){
                    if (first_base){
                        first_base = false;
                        out << name_line << "\n";
                    }
                    out << base;
                },
                [&](){ //non-ACGT base: finish outputting the line and create a new one
                    if (!first_base){
                        out << "\n";
                        first_base = true;
                        name_line = name_line + "^n";
                    }
                });
            if (!first_base){
                out << "\n";
            }
        });

        out.close();

        //append the chunk to the final output
        #pragma omp critical
        {
            std::ofstream final_out(output_file, std::ios::app | std::ios::binary);
            std::ifstream chunk_in(output_file_chunk, std::ios::binary);
            final_out << chunk_in.rdbuf();
            chunk_in.close();
            std::remove(output_file_chunk.c_str());
        }
    }
}

/**
 * @brief Go through the reads and, for every compressed kmer of the assembly, collect the uncompressed sequence(s) it
 * corresponds to. Candidates are counted per canonical compressed kmer, oriented on the canonical hash. A kmer is
 * "confirmed" as soon as one candidate is seen more than 3 times; otherwise the most frequent candidate is used
 * (rather than the first one seen, which carries the errors of whatever read happened to come first).
 * Both orientations are written at the very end, so that no thread/chunk can overwrite a confirmed sequence.
 */
// A "full" kmer is stored as "<o_0>,<o_1>,...,<o_{km-1}>\t<seq>", where o_j is the index in seq of the base of the
// j-th sampled position of the compressed kmer. This lets the expansion splice the ends of the contigs exactly at the
// sampled positions (no string matching, which fails in repeats or when the flank overhangs the rest of the contig).
static string encode_full(const vector<int>& offsets, const string& seq){
    string res;
    for (size_t j = 0 ; j < offsets.size() ; j++){
        if (j > 0) res += ",";
        res += std::to_string(offsets[j]);
    }
    return res + "\t" + seq;
}
static bool decode_full(const string& enc, vector<int>& offsets, string& seq){
    offsets.clear();
    size_t t = enc.find('\t');
    if (t == string::npos){ seq = enc; return false; }
    seq = enc.substr(t+1);
    std::stringstream ss(enc.substr(0, t));
    string o;
    while (std::getline(ss, o, ',')){
        offsets.push_back(std::stoi(o));
    }
    return true;
}
static string rc_full(const string& enc){
    if (enc == "") return "";
    vector<int> offsets; string seq;
    decode_full(enc, offsets, seq);
    vector<int> rc_offsets(offsets.size());
    for (size_t j = 0 ; j < offsets.size() ; j++){
        rc_offsets[j] = (int)seq.size() - 1 - offsets[offsets.size()-1-j];
    }
    return encode_full(rc_offsets, reverse_complement(seq));
}

// The candidate sequences are stored packed to save memory: 2 bits per base (or raw if the sequence contains a
// non-ACGT base), and the offsets of the full kmers as varint-encoded deltas
static void append_varint(string& out, uint32_t x){
    while (x >= 128){
        out += (char) (x % 128 + 128);
        x /= 128;
    }
    out += (char) x;
}
static uint32_t read_varint(const string& in, size_t& pos){
    uint32_t x = 0;
    uint32_t factor = 1;
    while ((unsigned char) in[pos] >= 128){
        x += ((unsigned char) in[pos] - 128) * factor;
        factor *= 128;
        pos++;
    }
    x += (unsigned char) in[pos] * factor;
    pos++;
    return x;
}
static void append_packed_sequence(string& out, const string& seq){
    bool only_ACGT = std::all_of(seq.begin(), seq.end(), is_ACGT);
    out += only_ACGT ? '\0' : '\1';
    append_varint(out, seq.size());
    if (!only_ACGT){
        out += seq;
        return;
    }
    unsigned char byte = 0;
    for (size_t i = 0 ; i < seq.size() ; i++){
        unsigned char code = seq[i] == 'A' ? 0 : (seq[i] == 'C' ? 1 : (seq[i] == 'G' ? 2 : 3));
        byte |= code << (2*(i%4));
        if (i%4 == 3){
            out += (char) byte;
            byte = 0;
        }
    }
    if (seq.size()%4 != 0){
        out += (char) byte;
    }
}
static string read_packed_sequence(const string& in, size_t& pos){
    bool only_ACGT = in[pos] == '\0';
    pos++;
    size_t length = read_varint(in, pos);
    if (!only_ACGT){
        pos += length;
        return in.substr(pos - length, length);
    }
    string seq (length, 'A');
    for (size_t i = 0 ; i < length ; i++){
        seq[i] = "ACGT"[((unsigned char) in[pos + i/4] >> (2*(i%4))) & 3];
    }
    pos += (length+3)/4;
    return seq;
}
static string pack_central(const string& seq){
    string packed;
    append_packed_sequence(packed, seq);
    return packed;
}
static string unpack_central(const string& packed){
    size_t pos = 0;
    return read_packed_sequence(packed, pos);
}
static string pack_full(const vector<int>& offsets, const string& seq){
    string packed;
    append_varint(packed, offsets.size());
    int previous = 0;
    for (int o : offsets){ //the offsets are increasing
        append_varint(packed, o - previous);
        previous = o;
    }
    append_packed_sequence(packed, seq);
    return packed;
}
//back to the "<o_0>,<o_1>,...\t<seq>" representation
static string unpack_full(const string& packed){
    size_t pos = 0;
    vector<int> offsets (read_varint(packed, pos));
    int previous = 0;
    for (int& o : offsets){
        o = previous + read_varint(packed, pos);
        previous = o;
    }
    return encode_full(offsets, read_packed_sequence(packed, pos));
}

struct Candidate {
    uint64_t hash;
    string packed;
    int count;
};

struct KmerCandidates {
    uint64_t other_hash = 0;          // non-canonical oriented hash
    bool need_central_canon = false, need_full_canon = false, need_central_other = false, need_full_other = false;
    bool confirmed_central = false, confirmed_full = false;
    string chosen_central, chosen_full; // packed, oriented on the canonical hash
    vector<Candidate> central_candidates, full_candidates;
    bool confirmed() const {
        return (confirmed_central || !(need_central_canon || need_central_other)) && (confirmed_full || !(need_full_canon || need_full_other));
    }
};

//the candidates are split in shards, each with its own lock, so that the threads do not all wait on the same lock
struct KmerShard {
    std::shared_mutex mtx;
    std::unordered_map<uint64_t, KmerCandidates> kmers; //canonical compressed hash -> candidates
};
static const int NUMBER_OF_KMER_SHARDS = 4096;

//register one observation of a (packed) candidate sequence; a candidate seen more than 3 times is confirmed
static void vote(vector<Candidate>& candidates, bool& confirmed, string& chosen, const string& packed){
    if (confirmed || packed == "") return;
    uint64_t hash = std::hash<string>()(packed);
    int count = 0;
    for (Candidate& c : candidates){
        if (c.hash == hash && c.packed == packed){
            count = ++c.count;
            break;
        }
    }
    if (count == 0){
        candidates.push_back({hash, packed, 1});
        count = 1;
    }
    if (count > 3){
        confirmed = true;
        chosen = packed;
        vector<Candidate>().swap(candidates); //free the memory
    }
}
//pick the most frequent candidate; ties are broken by taking the smallest sequence (unpacked), to be deterministic
static void choose(vector<Candidate>& candidates, bool confirmed, string& chosen, string (*unpack)(const string&)){
    if (confirmed) return;
    int best = 0;
    string best_unpacked;
    for (const Candidate& c : candidates){
        if (c.count > best){
            best = c.count;
            chosen = c.packed;
            best_unpacked = "";
        }
        else if (c.count == best){
            if (best_unpacked == ""){
                best_unpacked = unpack(chosen);
            }
            string unpacked = unpack(c.packed);
            if (unpacked < best_unpacked){
                chosen = c.packed;
                best_unpacked = unpacked;
            }
        }
    }
}

void go_through_the_reads_again_and_index_interesting_kmers(string reads_file,
    string assemblyFile,
    int order,
    int compression,
    int km,
    std::vector<uint64_t> &central_kmers_in_assembly,
    std::vector<uint64_t> &full_kmers_in_assembly,
    unordered_map<uint64_t, pair<unsigned long long,unsigned long long>>& kmers,
    string central_kmers_file,
    string full_kmers_file,
    int num_threads,
    bool homopolymer_compression){

    //one hash table lookup per oriented kmer instead of four binary searches
    const uint8_t CENTRAL = 1, FULL = 2;
    robin_hood::unordered_flat_map<uint64_t, uint8_t> needed_kmers;
    needed_kmers.reserve(central_kmers_in_assembly.size() + full_kmers_in_assembly.size());
    for (uint64_t h : central_kmers_in_assembly){
        needed_kmers[h] |= CENTRAL;
    }
    for (uint64_t h : full_kmers_in_assembly){
        needed_kmers[h] |= FULL;
    }
    auto flags_of = [&](uint64_t h) -> uint8_t {
        auto it = needed_kmers.find(h);
        return it == needed_kmers.end() ? 0 : it->second;
    };

    std::unique_ptr<KmerShard[]> shards (new KmerShard[NUMBER_OF_KMER_SHARDS]);
    auto shard_of = [&](uint64_t canonical_hash) -> KmerShard& {
        return shards[(canonical_hash >> 17) % NUMBER_OF_KMER_SHARDS];
    };

    ifstream input2(reads_file, std::ios::binary | std::ios::ate);
    std::streamoff file_size = input2.tellg();
    unsigned long long int size_of_chunk = 100000000;
    input2.close();

    std::atomic<unsigned long long> seq_num (0);
    std::atomic<unsigned long long> output_limit (0);
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int chunk = 0 ; chunk <= file_size/size_of_chunk ; chunk++){

        for_each_read_in_chunk(reads_file, chunk, size_of_chunk, [&](const string&, string& line){
            unsigned long long n = seq_num++;
            if (n >= output_limit && omp_get_thread_num() == 0){
                #pragma omp critical
                {
                    print_timestamp();
                    cout << " Processed " << n << " reads" << endl;
                }
                output_limit += 50000;
            }

            size_t pos_end_compressed = 0;
            long pos_begin_compressed = -km;
            long pos_middle_compressed = -6666;

            vector<int> positions_sampled (0);
            uint64_t hash_foward_compressed = 0, hash_reverse_compressed = 0;
            string compressed_read = "";
            compressed_read.reserve(line.size()/compression);

            sample_read(line, order, compression, homopolymer_compression, true,
                [&](char base, long pos_middle){
                    compressed_read += base;
                    roll(hash_foward_compressed, hash_reverse_compressed, km, compressed_read, pos_end_compressed, pos_begin_compressed, pos_middle_compressed, false);
                    positions_sampled.push_back(pos_middle);

                    if (positions_sampled.size() < km){
                        return;
                    }

                    uint8_t flags_fw = flags_of(hash_foward_compressed);
                    uint8_t flags_rv = flags_of(hash_reverse_compressed);
                    bool fc = flags_fw & CENTRAL, rc = flags_rv & CENTRAL;
                    bool ff = flags_fw & FULL, rf = flags_rv & FULL;
                    bool central_kmer_in_assembly = fc || rc;
                    bool full_kmer_in_assembly = ff || rf;

                    if (!central_kmer_in_assembly && !full_kmer_in_assembly){
                        return;
                    }

                    bool fw_is_canonical = hash_foward_compressed <= hash_reverse_compressed;
                    uint64_t canonical_hash = fw_is_canonical ? hash_foward_compressed : hash_reverse_compressed;
                    KmerShard& shard = shard_of(canonical_hash);

                    {
                        std::shared_lock<std::shared_mutex> lock(shard.mtx);
                        auto it = shard.kmers.find(canonical_hash);
                        if (it != shard.kmers.end() && it->second.confirmed()){
                            return;
                        }
                    }

                    //the sequences, oriented on the canonical hash and packed
                    string central_seq = "", full_seq = "";
                    if (central_kmer_in_assembly){
                        string seq = line.substr(positions_sampled[positions_sampled.size()-km+10], positions_sampled[positions_sampled.size()-1-10] - positions_sampled[positions_sampled.size()-km+10]+1);
                        central_seq = pack_central(fw_is_canonical ? seq : reverse_complement(seq));
                    }
                    if (full_kmer_in_assembly){
                        auto begin = std::max(0,positions_sampled[positions_sampled.size()-km] - order);
                        auto end = std::min((int)line.size()-1, positions_sampled[positions_sampled.size()-1] + order);
                        string seq = line.substr(begin, end - begin +1);
                        vector<int> offsets(km);
                        for (int j = 0 ; j < km ; j++){
                            offsets[j] = positions_sampled[positions_sampled.size()-km+j] - begin;
                        }
                        if (!fw_is_canonical){
                            vector<int> rc_offsets(km);
                            for (int j = 0 ; j < km ; j++){
                                rc_offsets[j] = (int)seq.size() - 1 - offsets[km-1-j];
                            }
                            offsets = rc_offsets;
                            seq = reverse_complement(seq);
                        }
                        full_seq = pack_full(offsets, seq);
                    }

                    std::unique_lock<std::shared_mutex> lock(shard.mtx);
                    KmerCandidates &kc = shard.kmers[canonical_hash];
                    kc.other_hash = fw_is_canonical ? hash_reverse_compressed : hash_foward_compressed;
                    kc.need_central_canon |= fw_is_canonical ? fc : rc;
                    kc.need_full_canon    |= fw_is_canonical ? ff : rf;
                    kc.need_central_other |= fw_is_canonical ? rc : fc;
                    kc.need_full_other    |= fw_is_canonical ? rf : ff;
                    vote(kc.central_candidates, kc.confirmed_central, kc.chosen_central, central_seq);
                    vote(kc.full_candidates, kc.confirmed_full, kc.chosen_full, full_seq);
                },
                [&](){ //non-ACGT base: start a new compressed read
                    positions_sampled.clear();
                    compressed_read.clear();
                    pos_end_compressed = 0;
                    pos_begin_compressed = -km;
                });
        });
    }

    //now choose a sequence for each kmer (most frequent candidate if not confirmed) and write both orientations
    std::ofstream out_central(central_kmers_file);
    std::ofstream out_full(full_kmers_file);
    for (int s = 0 ; s < NUMBER_OF_KMER_SHARDS ; s++){
        for (auto &p : shards[s].kmers){
            KmerCandidates &kc = p.second;
            choose(kc.central_candidates, kc.confirmed_central, kc.chosen_central, unpack_central);
            choose(kc.full_candidates, kc.confirmed_full, kc.chosen_full, unpack_full);
            string central_canon = kc.chosen_central == "" ? "" : unpack_central(kc.chosen_central);
            string full_canon = kc.chosen_full == "" ? "" : unpack_full(kc.chosen_full);
            string central_other = reverse_complement(central_canon), full_other = rc_full(full_canon);
            for (int orientation = 0 ; orientation < 2 ; orientation++){
                bool need_c = orientation == 0 ? kc.need_central_canon : kc.need_central_other;
                bool need_f = orientation == 0 ? kc.need_full_canon : kc.need_full_other;
                if (!need_c && !need_f) continue;
                string &c = orientation == 0 ? central_canon : central_other;
                string &f = orientation == 0 ? full_canon : full_other;
                unsigned long long position_central = 1, position_full = 1;
                if (need_c && c != ""){
                    position_central = out_central.tellp();
                    out_central << c << "\n";
                }
                if (need_f && f != ""){
                    position_full = out_full.tellp();
                    out_full << f << "\n";
                }
                kmers[orientation == 0 ? p.first : kc.other_hash] = {position_central, position_full};
            }
        }
        shards[s].kmers.clear();
    }
    out_central.close();
    out_full.close();
}


/**
 * @brief Takes in the reduced assembly and either list the kmers needed for expansion or actually expand the assembly
 *
 * @param mode either index or expand on wether you want to index the kmers needed for indexing or actually expand the assembly
 * @param asm_reduced reduced assembly
 * @param km size of k used for compression
 * @param compression compression factor (to estimate the length of the sequences that cannot be expanded)
 * @param central_kmers_needed list of the central kmers needed for expansion (initially empty for index mode, useless for expand mode)
 * @param full_kmers_needed list of the full kmers needed for expansion (initially empty for index mode, useless for expand mode)
 * @param central_kmers_file file containing central inflated kmers
 * @param full_kmers_file file containing full inflated kmers
 * @param kmers dictionary mapping compressed kmers to their positions in the kmer file
 * @param output output file (uselss for index mode)
 */
void expand_or_list_kmers_needed_for_expansion(string mode, string asm_reduced, int km, int compression, std::vector<uint64_t> &central_kmers_needed, std::vector<uint64_t> &full_kmers_needed, string central_kmers_file, string full_kmers_file, unordered_map<uint64_t, pair<unsigned long long,unsigned long long>>& kmers, string output){

    if (mode != "expand" && mode != "index"){
        cerr << "ERROR (code 190) mode not supported in expand_or_list_kmers_needed_for_expansion " << mode << "\n";
        exit(1);
    }

    ifstream input(asm_reduced);

    ofstream out;
    if (mode == "expand"){
        out.open(output);
        out << "H\tVN:Z:1.0\n";
    }

    //load the compressed sequences (the compressed assembly is small)
    unordered_map<std::string, std::string> compressed_sequences;
    string line;
    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            compressed_sequences[name] = sequence;
        }
    }
    input.close();

    //contig ends that are linked to another contig: there, the expanded sequence must stop exactly at the boundary sample
    std::unordered_set<std::string> linked_left, linked_right;
    //provide 10 bp left and right of all contigs if possible, based on links in the gfa, to improve expansion
    unordered_map<std::string, std::string> left_seq;
    unordered_map<std::string, std::string> right_seq;
    input.open(asm_reduced);
    while (std::getline(input, line))
    {
        if (line[0] == 'L')
        {
            string contig1;
            string contig2;
            string orientation1;
            string orientation2;
            string dont_care;
            string cigar;
            std::stringstream ss(line);
            ss >> dont_care >> contig1 >> orientation1 >> contig2 >> orientation2 >> cigar;

            //parse the cigar and output an error if there is something else than M
            int length_of_overlap;
            if (cigar.find("M") == std::string::npos){
                cerr << "ERROR (code 322) cigar not supported " << cigar << "\n";
                exit(1);
            }
            else{
                length_of_overlap = std::stoi(cigar.substr(0, cigar.find("M")));
            }

            auto it1 = compressed_sequences.find(contig1);
            auto it2 = compressed_sequences.find(contig2);
            if (it1 == compressed_sequences.end() || it2 == compressed_sequences.end()){
                cerr << "WARNING (code 323) link between unknown contigs ignored: " << line << "\n";
                continue;
            }

            if (orientation1 == "+") linked_right.insert(contig1); else linked_left.insert(contig1);
            if (orientation2 == "+") linked_left.insert(contig2); else linked_right.insert(contig2);

            const string& sequence1 = it1->second;
            string first_10_seq1, last_10_seq1;
            if (sequence1.size() > 10+length_of_overlap){
                first_10_seq1 = sequence1.substr(length_of_overlap, 10);
                last_10_seq1 = sequence1.substr(sequence1.size()-10-length_of_overlap, 10);
            }

            const string& sequence2 = it2->second;
            string first_10_seq2, last_10_seq2;
            if (sequence2.size() > 10+length_of_overlap){
                first_10_seq2 = sequence2.substr(length_of_overlap, 10);
                last_10_seq2 = sequence2.substr(sequence2.size()-10-length_of_overlap, 10);
            }

            if (orientation1 == "+"){
                if (orientation2 == "+"){
                    if (first_10_seq2 != ""){
                        right_seq[contig1] = first_10_seq2;
                    }
                    if (last_10_seq1 != ""){
                        left_seq[contig2] = last_10_seq1;
                    }
                }
                else{
                    if (last_10_seq2 != ""){
                        right_seq[contig1] = reverse_complement(last_10_seq2);
                    }
                    if (last_10_seq1 != ""){
                        right_seq[contig2] = reverse_complement(last_10_seq1);
                    }
                }
            }
            else{
                if (orientation2 == "+"){
                    if (first_10_seq2 != ""){
                        left_seq[contig1] = reverse_complement(first_10_seq2);
                    }
                    if (first_10_seq1 != ""){
                        left_seq[contig2] = reverse_complement(first_10_seq1);
                    }
                }
                else{
                    if (last_10_seq2 != ""){
                        left_seq[contig1] = last_10_seq2;
                    }
                    if (first_10_seq1 != ""){
                        right_seq[contig2] = first_10_seq1;
                    }
                }
            }
        }
    }
    input.close();

    //open the kmers file
    ifstream central_kmers_input(central_kmers_file); //not used in index mode
    ifstream full_kmers_input(full_kmers_file); //not used in index mode

    input.open(asm_reduced);
    int number_of_missing_kmers = 0;
    string line2;

    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;

            //extend left and right if possible
            auto left_it = left_seq.find(name);
            auto right_it = right_seq.find(name);
            bool has_left_ext = left_it != left_seq.end();
            bool has_right_ext = right_it != right_seq.end();
            if (has_left_ext){
                sequence = left_it->second + sequence;
            }
            if (has_right_ext){
                sequence += right_it->second;
            }

            if (sequence.size() < km){ //too short to be expanded
                if (mode == "expand"){
                    cerr << "WARNING (code 743) sequence too short \n";
                    string expanded_sequence (sequence.size() * compression, 'N');
                    out << "S\t" << name << "\t" << expanded_sequence;
                    while (ss >> sequence){
                        out << "\t" << sequence;
                    }
                    out << "\n";
                }
                continue;
            }

            //expand the sequence (focusing on the central part of each kmer and thus missing the two ends)
            //central part of the kmer starting at i = samples i+10 .. i+km-11 (both included)
            string expanded_sequence = "";
            int i = 0;
            int last_i = 0;
            int length_of_central_kmers = 10*compression+1; //central part = 11 samples, i.e. ~10*compression bases (refined with the actual lengths)
            for (i = 0; i <= (int)sequence.size()-km; i+= km-20-1){ //-20 because we only take the central part of each kmer
                last_i = i;
                string kmer = sequence.substr(i, km);
                uint64_t hash_foward_kmer = hash_string(km, kmer, false);
                if (mode == "index"){
                    central_kmers_needed.push_back(hash_foward_kmer);
                }
                else{
                    auto it = kmers.find(hash_foward_kmer);
                    if (it != kmers.end() && it->second.first != 1){
                        central_kmers_input.seekg(it->second.first);
                        std::getline(central_kmers_input, line2);
                        const string& central_seq = line2;

                        if (central_seq.size() > 1){
                            length_of_central_kmers = central_seq.size();
                        }
                        if (i == 0 || central_seq.size() == 0){
                            expanded_sequence += central_seq;
                        }
                        else{ //don't append the first base, it has already been appended in the previous kmer
                            expanded_sequence.append(central_seq.begin()+1, central_seq.end());
                        }
                    }
                    else{
                        number_of_missing_kmers++;
                        //same length as an expanded central part (whose first base is shared with the previous one)
                        expanded_sequence += string(i == 0 ? length_of_central_kmers : length_of_central_kmers-1, 'N');
                    }
                }
            }

            // create the beginning of the sequence if there was no left extension: the expanded sequence currently
            // starts at the base of sample 10 of the first kmer; prepend what is before, using the sample offsets
            if (!has_left_ext){
                string first_kmer = sequence.substr(0, km);
                uint64_t hash_foward_kmer = hash_string(km, first_kmer, false);
                if (mode == "index"){
                    full_kmers_needed.push_back(hash_foward_kmer);
                }
                else {
                    vector<int> off; string full_kmer;
                    auto it = kmers.find(hash_foward_kmer);
                    if (it != kmers.end() && it->second.second != 1){
                        full_kmers_input.seekg(it->second.second);
                        std::getline(full_kmers_input, line2);
                        decode_full(line2, off, full_kmer);
                    }
                    if ((int)off.size() == km){
                        //if the contig is linked on its left, start exactly at the first sampled base (no overhang)
                        int cut = linked_left.count(name) ? off[0] : 0;
                        expanded_sequence = full_kmer.substr(cut, off[10]-cut) + expanded_sequence;
                    }
                    else{
                        number_of_missing_kmers++;
                    }
                }
            }

            // finish the sequence: the expanded sequence currently ends at the base of sample (last_i+km-11) of the
            // (extended) sequence; append the rest using the sample offsets of the last full kmer
            int j = (int)sequence.size()-km; //start of the last kmer
            string last_kmer = sequence.substr(j, km);
            uint64_t hash_foward_kmer = hash_string(km, last_kmer, false);

            if (mode == "index"){
                full_kmers_needed.push_back(hash_foward_kmer);
            }
            else
            {
                vector<int> off; string full_kmer;
                auto it = kmers.find(hash_foward_kmer);
                if (it != kmers.end() && it->second.second != 1){
                    full_kmers_input.seekg(it->second.second);
                    std::getline(full_kmers_input, line2);
                    decode_full(line2, off, full_kmer);
                }
                if ((int)off.size() == km){
                    int e = last_i + km - 11 - j; //sample of the last kmer at which expanded_sequence currently ends (inclusive)
                    int stop; //position in full_kmer after the last base to keep
                    if (has_right_ext){
                        stop = off[km-11]+1; //the last 10 samples belong to the next contig
                    }
                    else if (linked_right.count(name)){
                        stop = off[km-1]+1; //do not overhang beyond the last sampled base
                    }
                    else{
                        stop = full_kmer.size(); //tip: keep the whole flank
                    }
                    if (stop > off[e]+1){
                        expanded_sequence += full_kmer.substr(off[e]+1, stop-off[e]-1);
                    }
                }
                else{
                    number_of_missing_kmers++;
                }

                out << "S\t" << name << "\t" << expanded_sequence;
                while (ss >> sequence){
                    out << "\t" << sequence;
                }
                out << "\n";
            }
        }
        else if (mode == "expand" && line[0] != 'H'){
            out << line << "\n";
        }
    }
    if (mode == "expand"){
            cout << "WARNING for developers: Number of missing kmers " << number_of_missing_kmers << endl;
            out.close();
    }
    central_kmers_input.close();
    full_kmers_input.close();

}
