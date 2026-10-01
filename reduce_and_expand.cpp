#include "reduce_and_expand.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <unordered_set>
#include <chrono>
#include <thread>
#include <algorithm>
#include <omp.h> //for efficient parallelization
#include <filesystem>
// #include <zlib.h>  // Include the gzstream header
#include <set>
#include <map>
#include <shared_mutex>
#include <mutex>


#include "robin_hood.h"
#include "basic_graph_manipulation.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using robin_hood::unordered_map;
using std::pair;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::unordered_set;

/**
 * @brief MSR the input sequencing
 * 
 * @param input_file 
 * @param output_file 
 * @param sample_input_file Subsample of the input reads that strive to represent all sequences
 * @param context_length 
 * @param compression 
 * @param km 
 * @param min_abundance 
 * @param kmers maps a kmer to the uncompressed seq
 **/
void reduce(string input_file, string output_file, int order, int compression, int num_threads, bool homopolymer_compression) {

    time_t now2 = time(0);
    tm *ltm2 = localtime(&now2);
    cout << "[" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "] Starting pipeline" << endl;


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

    unsigned long long int seq_num = 0;
    long long output_limit = 0;
    //parallelize on num_threads threads
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int chunk = 0 ; chunk <= file_size/size_of_chunk ; chunk++){

        std::ifstream input(input_file);
        input.seekg(chunk*size_of_chunk);
    
        string output_file_chunk = output_file + "_"+ std::to_string(chunk);
        std::ofstream out(output_file_chunk);
        std::string line;
        bool next_line_is_seq = false;
        string name_line = "";
        unordered_map<uint64_t, short> number_of_kmer_occurences; //count the number of occurence of some kmers to know if the read is in an "already seen" context

        while (std::getline(input, line))
        {
            if (line[0] == '>')
            {
                //out << line << "\n";
                if (seq_num > output_limit && omp_get_thread_num() == 0){
                    #pragma omp critical
                    {
                        //display the date and time
                        time_t now = time(0);
                        tm *ltm = localtime(&now);
                        cout << "[" << 1+ ltm->tm_mday << "/" << 1 + ltm->tm_mon << "/" << 1900 + ltm->tm_year << " " << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec << "]" << " Compressed " << seq_num << " reads" << endl;
                        output_limit += 50000;
                    }
                }
                seq_num++;
                next_line_is_seq = true;
                name_line = line;
            }
            else if (next_line_is_seq){
                //let's launch the foward and reverse rolling hash
                uint64_t hash_foward = 0;
                uint64_t hash_reverse = 0;
                size_t pos_end = 0;
                long pos_middle = -6666; //special value to not compute the middle position
                long pos_begin = -order;
                bool first_base = true; //to output the name line when the first base is outputted (not before to avoid empty lines)
                next_line_is_seq = false;
                int number_of_hashed_bases = 0;
                bool this_sequence_has_not_been_seen_before = false; //to keep only "new" sequences to accelerate decompression in the future

                while (roll(hash_foward, hash_reverse, order, line, pos_end, pos_begin, pos_middle, homopolymer_compression)){
                    if (line[pos_end] != 'A' && line[pos_end] != 'C' && line[pos_end] != 'G' && line[pos_end] != 'T'){
                        //then finish outputting the line and create a new one
                        if (!first_base){
                            out << "\n";
                            first_base = true;
                            name_line = name_line + "^n";
                        }
                        //reinitialize the rolling hash
                        number_of_hashed_bases = 0;

                    }
                    if (number_of_hashed_bases >= order){
                        if (hash_foward<hash_reverse && hash_foward % compression == 0){
                            if (first_base){
                                first_base = false;
                                out << name_line << "\n";
                            }
                            out << "ACGT"[(hash_foward/compression)%4];

                            if ((hash_foward / compression) % 10 == 0){
                                // if (name_line == ">SRR13128013.1 1 length=2993"){
                                //     cout << "hash fw in read SRR13128013.1 1 length=2993 " << hash_foward << " for kmer " << line.substr(pos_begin, pos_end-pos_begin) << "\n";
                                // }
                                // else if (name_line == ">SRR13128013.26581 26581 length=7693"){
                                //     cout << "hash fw in read SRR13128013.20422 20422 length=16769 " << hash_foward << "\n";
                                // }
                                #pragma omp critical
                                {
                                    number_of_kmer_occurences[hash_foward]++;
                                    if (number_of_kmer_occurences[hash_foward] == 2){ //if it is exactly the second time we see this kmer, then we know that this sequence has not been seen before (once is error)
                                        this_sequence_has_not_been_seen_before = true;
                                    }
                                    else if (this_sequence_has_not_been_seen_before == true){
                                        number_of_kmer_occurences[hash_foward]++; //to make sure it is at least 2, i.e. does not need to be kept in another read
                                    }
                                }
                            }
                        }
                        else if (hash_foward>=hash_reverse && hash_reverse % compression == 0){
                            if (first_base){
                                first_base = false;
                                out << name_line << "\n";
                            }
                            out << "TGCA"[(hash_reverse/compression)%4];

                            // if ((hash_reverse / compression) % 10 == 0){
                            //     // if (name_line == ">SRR13128013.1 1 length=2993"){
                            //     //     cout << "hash rv in read SRR13128013.1 1 length=2993 " << hash_reverse << "\n";
                            //     // }
                            //     // else if (name_line == ">SRR13128013.26581 26581 length=7693"){
                            //     //     cout << "hash rv in read SRR13128013.20422 20422 length=16769 " << hash_reverse << "\n";
                            //     // }
                            //     #pragma omp critical
                            //     {
                            //         number_of_kmer_occurences[hash_reverse]++;
                            //         if (number_of_kmer_occurences[hash_reverse] == 2){ //if it is exactly the second time we see this kmer, then we know that this sequence has not been seen before (once is error)
                            //             this_sequence_has_not_been_seen_before = true;
                            //         }
                            //         else if (this_sequence_has_not_been_seen_before == true){
                            //             number_of_kmer_occurences[hash_reverse]++; //to make sure it is at least 2, i.e. does not need to be kept in another read
                            //         }
                            //     }
                            // }
                        }
                    }
                    number_of_hashed_bases++;
                }
                out << "\n";

                if (input.tellg() > (chunk+1)*size_of_chunk){
                    break;
                }
            }
        }

        input.close();
        out.close();

        //append the chunk to the final output
        #pragma omp critical
        {
            system(("cat " + output_file_chunk + " >> " + output_file).c_str());
            system(("rm " + output_file_chunk).c_str());
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

struct KmerCandidates {
    uint64_t other_hash = 0;          // non-canonical oriented hash
    bool need_central_canon = false, need_full_canon = false, need_central_other = false, need_full_other = false;
    bool confirmed_central = false, confirmed_full = false;
    string chosen_central, chosen_full; // oriented on the canonical hash
    std::map<string, int> central_candidates, full_candidates;
    bool confirmed() const {
        return (confirmed_central || !(need_central_canon || need_central_other)) && (confirmed_full || !(need_full_canon || need_full_other));
    }
};

//register one observation of a candidate sequence; a candidate seen more than 3 times is confirmed
static void vote(std::map<string,int>& candidates, bool& confirmed, string& chosen, const string& seq){
    if (confirmed || seq == "") return;
    int count = ++candidates[seq];
    if (count > 3){
        confirmed = true;
        chosen = seq;
        candidates.clear();
    }
}
//pick the most frequent candidate (deterministic tie-breaking thanks to std::map)
static void choose(std::map<string,int>& candidates, bool confirmed, string& chosen){
    if (confirmed) return;
    int best = 0;
    for (auto &c : candidates){
        if (c.second > best){
            best = c.second;
            chosen = c.first;
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

    std::unordered_map<uint64_t, KmerCandidates> all_kmers; //canonical compressed hash -> candidates
    std::shared_mutex mtx;

    ifstream input2(reads_file, std::ios::binary | std::ios::ate);
    std::streamoff file_size = input2.tellg();
    unsigned long long int size_of_chunk = 100000000;
    input2.close();

    int seq_num = 0;
    long long output_limit = 0;
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int chunk = 0 ; chunk <= file_size/size_of_chunk ; chunk++){

        std::ifstream input(reads_file);
        input.seekg(chunk*size_of_chunk);
    
        std::string line;
        bool next_line_is_seq = false;

        while (std::getline(input, line)){

            if (line[0] == '>')
            {
                if (seq_num > output_limit && omp_get_thread_num() == 0){
                    #pragma omp critical
                    {
                        time_t now = time(0);
                        tm *ltm = localtime(&now);
                        cout << "[" << 1 + ltm->tm_mday << "/" << 1 + ltm->tm_mon << "/" << 1900 + ltm->tm_year << " " << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec << "]" << " Processed " << seq_num << " reads" << endl;
                        output_limit += 50000;
                    }
                }
                #pragma omp atomic
                seq_num++;
                next_line_is_seq = true;
            }
            else if (next_line_is_seq){
                uint64_t hash_foward = 0;
                uint64_t hash_reverse = 0;
                size_t pos_end = 0;
                long pos_begin = -order;
                long pos_middle = -(order+1)/2;
                next_line_is_seq = false;
                int number_of_hashed_bases = 0;

                size_t pos_end_compressed = 0;
                long pos_begin_compressed = -km;
                long pos_middle_compressed = -6666;

                vector<int> positions_sampled (0);
                uint64_t hash_foward_compressed = 0, hash_reverse_compressed = 0;
                string compressed_read = "";
                compressed_read.reserve(line.size()/compression);

                while (roll(hash_foward, hash_reverse, order, line, pos_end, pos_begin, pos_middle, homopolymer_compression)){
                    if (number_of_hashed_bases >= order){

                        if (line[pos_end] != 'A' && line[pos_end] != 'C' && line[pos_end] != 'G' && line[pos_end] != 'T'){
                            positions_sampled.clear();
                            number_of_hashed_bases = 0;
                        }

                        if ((hash_foward<hash_reverse && hash_foward % compression == 0) || (hash_foward>=hash_reverse && hash_reverse % compression == 0)){

                            if (hash_foward<hash_reverse){
                                compressed_read += "ACGT"[(hash_foward/compression)%4];
                            }
                            else{
                                compressed_read += "TGCA"[(hash_reverse/compression)%4];
                            }

                            roll(hash_foward_compressed, hash_reverse_compressed, km, compressed_read, pos_end_compressed, pos_begin_compressed, pos_middle_compressed, false);
                            positions_sampled.push_back(pos_middle);

                            if (positions_sampled.size() >= km){

                                bool fc = std::binary_search(central_kmers_in_assembly.begin(), central_kmers_in_assembly.end(), hash_foward_compressed);
                                bool rc = std::binary_search(central_kmers_in_assembly.begin(), central_kmers_in_assembly.end(), hash_reverse_compressed);
                                bool ff = std::binary_search(full_kmers_in_assembly.begin(), full_kmers_in_assembly.end(), hash_foward_compressed);
                                bool rf = std::binary_search(full_kmers_in_assembly.begin(), full_kmers_in_assembly.end(), hash_reverse_compressed);
                                bool central_kmer_in_assembly = fc || rc;
                                bool full_kmer_in_assembly = ff || rf;

                                if (central_kmer_in_assembly || full_kmer_in_assembly){
                                
                                    bool fw_is_canonical = hash_foward_compressed <= hash_reverse_compressed;
                                    uint64_t canonical_hash = fw_is_canonical ? hash_foward_compressed : hash_reverse_compressed;

                                    bool already_confirmed = false;
                                    {
                                        std::shared_lock<std::shared_mutex> lock(mtx);
                                        auto it = all_kmers.find(canonical_hash);
                                        already_confirmed = (it != all_kmers.end() && it->second.confirmed());
                                    }
                                    if (!already_confirmed){

                                    string central_seq = "", full_seq = "";
                                    if (central_kmer_in_assembly){
                                        central_seq = line.substr(positions_sampled[positions_sampled.size()-km+10], positions_sampled[positions_sampled.size()-1-10] - positions_sampled[positions_sampled.size()-km+10]+1);
                                    }
                                    if (full_kmer_in_assembly){
                                        auto begin = std::max(0,positions_sampled[positions_sampled.size()-km] - order);
                                        auto end = std::min((int)line.size()-1, positions_sampled[positions_sampled.size()-1] + order);
                                        vector<int> offsets(km);
                                        for (int j = 0 ; j < km ; j++){
                                            offsets[j] = positions_sampled[positions_sampled.size()-km+j] - begin;
                                        }
                                        full_seq = encode_full(offsets, line.substr(begin, end - begin +1));
                                    }
                                    //orient the sequences on the canonical hash
                                    if (!fw_is_canonical){
                                        central_seq = reverse_complement(central_seq);
                                        full_seq = rc_full(full_seq);
                                    }

                                    {
                                        std::unique_lock<std::shared_mutex> lock(mtx);
                                        KmerCandidates &kc = all_kmers[canonical_hash];
                                        kc.other_hash = fw_is_canonical ? hash_reverse_compressed : hash_foward_compressed;
                                        kc.need_central_canon |= fw_is_canonical ? fc : rc;
                                        kc.need_full_canon    |= fw_is_canonical ? ff : rf;
                                        kc.need_central_other |= fw_is_canonical ? rc : fc;
                                        kc.need_full_other    |= fw_is_canonical ? rf : ff;
                                        vote(kc.central_candidates, kc.confirmed_central, kc.chosen_central, central_seq);
                                        vote(kc.full_candidates, kc.confirmed_full, kc.chosen_full, full_seq);
                                    }
                                    } //end if !already_confirmed
                                }
                            }
                        }
                    }
                    number_of_hashed_bases++;
                }
                if (input.tellg() > (chunk+1)*size_of_chunk){
                    break;
                }
            } 
        }
        input.close();
    }

    //now choose a sequence for each kmer (most frequent candidate if not confirmed) and write both orientations
    std::ofstream out_central(central_kmers_file);
    std::ofstream out_full(full_kmers_file);
    for (auto &p : all_kmers){
        KmerCandidates &kc = p.second;
        choose(kc.central_candidates, kc.confirmed_central, kc.chosen_central);
        choose(kc.full_candidates, kc.confirmed_full, kc.chosen_full);
        string central_canon = kc.chosen_central, full_canon = kc.chosen_full;
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
    out_central.close();
    out_full.close();
}


/**
 * @brief Takes in the reduced assembly and either list the kmers needed for expansion or actually expand the assembly
 * 
 * @param mode either index or expand on wether you want to index the kmers needed for indexing or actually expand the assembly
 * @param asm_reduced reduced assembly
 * @param km size of k used for compression
 * @param central_kmers_needed list of the central kmers needed for expansion (initially empty for index mode, useless for expand mode)
 * @param full_kmers_needed list of the full kmers needed for expansion (initially empty for index mode, useless for expand mode)
 * @param central_kmers_file file containing central inflated kmers
 * @param full_kmers_file file containing full inflated kmers
 * @param kmers dictionary mapping compressed kmers to their positions in the kmer file
 * @param output output file (uselss for index mode)
 */
void expand_or_list_kmers_needed_for_expansion(string mode, string asm_reduced, int km, std::vector<uint64_t> &central_kmers_needed, std::vector<uint64_t> &full_kmers_needed, string central_kmers_file, string full_kmers_file, unordered_map<uint64_t, pair<unsigned long long,unsigned long long>>& kmers, string output){
    
    if (mode != "expand" && mode != "index"){
        cerr << "ERROR (code 190) mode not supported in expand_or_list_kmers_needed_for_expansion " << mode << "\n";
        exit(1);
    }
    
    ifstream input(asm_reduced);

    ofstream out(output);
    if (mode == "expand"){
        ofstream out(output);
        out << "H\tVN:Z:1.0\n";
    }

    //first index the first and last 10 bases of all contigs
    unordered_map<std::string, std::string> first_10;
    unordered_map<std::string, std::string> last_10;
    unordered_map<std::string, long long> contig_locations;
    long long current_position = 0;
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
            contig_locations[name] = current_position;
        }
        current_position += line.size() + 1;
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

            if (orientation1 == "+") linked_right.insert(contig1); else linked_left.insert(contig1);
            if (orientation2 == "+") linked_left.insert(contig2); else linked_right.insert(contig2);

            //record the position in the file to come back after having retrieved the sequences
            long long position_in_file = input.tellg();

            //retrieve the sequences of the contigs
            input.seekg(contig_locations[contig1]);
            std::getline(input, line);
            string sequence1;
            std::stringstream ss1(line);
            ss1 >> dont_care >> contig1 >> sequence1;
            string first_10_seq1, last_10_seq1;
            if (sequence1.size() > 10+length_of_overlap){
                first_10_seq1 = sequence1.substr(length_of_overlap, 10);
                last_10_seq1 = sequence1.substr(sequence1.size()-10-length_of_overlap, 10);
            }

            input.seekg(contig_locations[contig2]);
            std::getline(input, line);
            string sequence2;
            std::stringstream ss2(line);
            ss2 >> dont_care >> contig2 >> sequence2;
            string first_10_seq2, last_10_seq2;
            if (sequence2.size() > 10+length_of_overlap){
                first_10_seq2 = sequence2.substr(length_of_overlap, 10);
                last_10_seq2 = sequence2.substr(sequence2.size()-10-length_of_overlap, 10);
            }

            //go back to the position in the file
            input.seekg(position_in_file);

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
    long long position_in_file;
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
            if (left_seq.find(name) != left_seq.end()){
                sequence = left_seq[name] + sequence;
            }
            if (right_seq.find(name) != right_seq.end()){
                sequence += right_seq[name];
            }

            if (sequence.size() < km && mode == "expand"){
                cerr << "WARNING (code 743) sequence too short \n";
                string expanded_sequence (sequence.size() * 20, 'N');
                out << "S\t" << name << "\t" << expanded_sequence;
                while (ss >> sequence){
                    out << "\t" << sequence;
                }
                out << "\n";
                continue;                
            }
            
            //expand the sequence (focusing on the central part of each kmer and thus missing the two ends)
            //central part of the kmer starting at i = samples i+10 .. i+km-11 (both included)
            string expanded_sequence = "";
            int i = 0;
            int last_i = 0;
            int length_of_central_kmers = 1 ;
            for (i = 0; i <= (int)sequence.size()-km; i+= km-20-1){ //-20 because we only take the central part of each kmer
                last_i = i;
                string kmer = sequence.substr(i, km);
                uint64_t hash_foward_kmer = hash_string(km, kmer, false);
                if (mode == "index"){
                    central_kmers_needed.push_back(hash_foward_kmer);
                }
                else{
                    if (kmers.find(hash_foward_kmer) != kmers.end() && kmers[hash_foward_kmer].first != 1){
                        central_kmers_input.seekg(kmers[hash_foward_kmer].first);
                        std::getline(central_kmers_input, line2);
                        string central_seq = line2;

                        length_of_central_kmers = central_seq.size();
                        if (i == 0 || central_seq.size() == 0){
                            expanded_sequence += central_seq;
                        }
                        else{ //don't append the first base, it has already been appended in the previous kmer
                            expanded_sequence.append(central_seq.begin()+1, central_seq.end());
                        }
                    }
                    else{
                        number_of_missing_kmers++;
                        expanded_sequence += string(length_of_central_kmers, 'N');
                    }
                }
            }

            bool has_left_ext = left_seq.find(name) != left_seq.end();
            bool has_right_ext = right_seq.find(name) != right_seq.end();

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
                    if (kmers.find(hash_foward_kmer) != kmers.end() && kmers[hash_foward_kmer].second != 1){
                        full_kmers_input.seekg(kmers[hash_foward_kmer].second);
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
                if (kmers.find(hash_foward_kmer) != kmers.end() && kmers[hash_foward_kmer].second != 1){
                    full_kmers_input.seekg(kmers[hash_foward_kmer].second);
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
        else if (mode == "expand"){
            out << line << "\n";
        }
    }
    if (mode == "expand"){
            cout << "WARNING for developers: Number of missing kmers " << number_of_missing_kmers << endl;
    }
    else{
        // cout << "Number of central kmers needed " << central_kmers_needed.size() << endl;
        // cout << "Number of full kmers needed " <<  full_kmers_needed.size() << endl;
    }
    out.close();
    central_kmers_input.close();
    full_kmers_input.close();

}
















