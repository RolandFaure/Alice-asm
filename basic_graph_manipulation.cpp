#include "basic_graph_manipulation.h"
#include "robin_hood.h"

#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <chrono>
#include <omp.h>
#include <memory>
#include <cmath>
#include <atomic>
#include <functional>

using std::cout;
using std::endl;
using std::string;
using std::set;
using robin_hood::unordered_flat_map;
using robin_hood::unordered_map;
using std::cerr;
using std::pair;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::unordered_set;
using std::make_pair;
using std::max;
using std::min;
using std::stringstream;

string reverse_complement(string& seq){
    string rc (seq.size(), 'N');
    for (int i = seq.size() - 1 ; i >= 0; i--){
        switch (seq[i]){
            case 'A':
                rc[seq.size() - 1 - i] = 'T';
                break;
            case 'T':
                rc[seq.size() - 1 - i] = 'A';
                break;
            case 'C':
                rc[seq.size() - 1 - i] = 'G';
                break;
            case 'G':
                rc[seq.size() - 1 - i] = 'C';
                break;
            default:
                rc[seq.size() - 1 - i] = 'N';
                break;
        }
    }
    return rc;
}

/**
 * @brief Read the end of a contig (as it is oriented in a link) directly in the gfa file, without loading the whole contig
 * 
 * @param gfa file opened on the gfa
 * @param seq_location position in the file of the first base of the sequence of the contig, and length of the sequence
 * @param orientation orientation of the contig in the link
 * @param suffix true to get the last bases of the oriented contig, false to get the first bases
 * @param length number of bases wanted (less if the contig is shorter)
 */
static string read_oriented_end(ifstream& gfa, const pair<long int, long int>& seq_location, const string& orientation, bool suffix, long int length){
    length = min(length, seq_location.second);
    //the suffix of the reverse complement is the reverse complement of the prefix
    bool read_suffix = (suffix == (orientation == "+"));
    long int start = read_suffix ? seq_location.first + seq_location.second - length : seq_location.first;
    string seq (length, 'N');
    gfa.clear();
    gfa.seekg(start);
    gfa.read(&seq[0], length);
    if (orientation == "-"){
        seq = reverse_complement(seq);
    }
    return seq;
}

/**
 * @brief Number of bases of the second segment covered by the overlap described by a GFA CIGAR ("*" means 0)
 */
int overlap_length_of_CIGAR(const string& cigar){
    int length = 0;
    int number = 0;
    for (char c : cigar){
        if (c >= '0' && c <= '9'){
            number = 10*number + (c - '0');
        }
        else {
            if (c == 'M' || c == '=' || c == 'X' || c == 'I'){
                length += number;
            }
            number = 0;
        }
    }
    return length;
}

long long total_sequence_length(std::string gfa){
    ifstream input(gfa);
    string line;
    long long total = 0;
    while (std::getline(input, line)){
        if (line[0] == 'S'){
            size_t seq_start = line.find('\t', line.find('\t') + 1) + 1;
            total += std::min(line.find('\t', seq_start), line.size()) - seq_start;
        }
    }
    return total;
}

/**
 * @brief Recompute the overlaps of the expanded graph, whose L lines still carry the overlaps of the compressed graph
 * 
 * @param max_overlap maximum overlap searched
 * @param default_overlap overlap used if no overlap is found
 * @param bases_per_compressed_base average number of bases per compressed base, to estimate the expected overlap
 */
void compute_exact_CIGARs(std::string gfa_in, std::string gfa_out, int max_overlap, int default_overlap, double bases_per_compressed_base, int num_threads){

    //go through the graph and for all links, compute the exact CIGAR (that will be only M)
    ifstream input(gfa_in);
    ofstream out(gfa_out);

    string line;
    //first index the position and length of the sequence of every contig in the file
    unordered_map<string, pair<long int, long int>> seq_location;
    long int pos = 0;
    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            size_t name_start = line.find('\t') + 1;
            size_t seq_start = line.find('\t', name_start) + 1;
            size_t seq_end = std::min(line.find('\t', seq_start), line.size());
            string name = line.substr(name_start, seq_start - 1 - name_start);
            seq_location[name] = {pos + (long int) seq_start, (long int) (seq_end - seq_start)};
        }
        pos += line.size() + 1;
    }

    input.close();
    input.open(gfa_in);

    //each thread reads the contig ends with its own stream
    vector<std::unique_ptr<ifstream>> sequence_readers;
    for (int t = 0 ; t < num_threads ; t++){
        sequence_readers.emplace_back(new ifstream(gfa_in, std::ios::binary));
    }

    //compute the overlap of one L line (returns the line to output, or "" to drop the link)
    auto compute_link = [&](const string& link_line, ifstream& reader) -> string {
        string name1;
        string name2;
        string orientation1;
        string orientation2;
        string dont_care;
        string compressed_cigar;
        std::stringstream ss(link_line);
        ss >> dont_care >> name1 >> orientation1 >> name2 >> orientation2 >> compressed_cigar;

        auto it1 = seq_location.find(name1);
        auto it2 = seq_location.find(name2);
        if (it1 == seq_location.end() || it2 == seq_location.end()){
            #pragma omp critical
            cerr << "WARNING: link between unknown contigs ignored: " << link_line << "\n";
            return "";
        }

        //the overlap spans the (compressed overlap) shared samples: estimate its length in bases
        int compressed_overlap = overlap_length_of_CIGAR(compressed_cigar);
        double expected_overlap = compressed_overlap > 0 ? (compressed_overlap-1)*bases_per_compressed_base + 1 : default_overlap;

        //get the end of the first contig and the beginning of the second one (only the max_overlap bases that can overlap)
        const auto& location1 = it1->second;
        const auto& location2 = it2->second;
        long int size1 = location1.second;
        long int size2 = location2.second;
        string seq1 = read_oriented_end(reader, location1, orientation1, true, max_overlap);
        string seq2 = read_oriented_end(reader, location2, orientation2, false, max_overlap);

        int overlap = 30;
        if (size1 < 30 || size2 < 30){ //can happen if they could not be reconstructed
            overlap = 0;
        }
        else{
            //among the overlaps >= 30 for which the whole suffix of seq1 equals the prefix of seq2, take the one
            //closest to the expected overlap (in repeats / low-complexity sequence there can be several, e.g. one per
            //period of a tandem repeat). If no exact overlap exists (e.g. an error at the end of a contig), fall back
            //on the overlap closest to the expected one for which the 30bp seed matches
            string end1 = seq1.substr(seq1.size()-30, 30); 
            long int max_possible = min((long int) max_overlap, min(size1, size2));
            auto closer = [&](int candidate, int current){
                return current == -1 || std::abs(candidate - expected_overlap) < std::abs(current - expected_overlap);
            };
            //scan the overlaps from `from` to `to` (included)
            auto scan = [&](long int from, long int to, int& best_exact, int& best_seed_hit){
                for (long int o = std::max(30L, from) ; o <= std::min(to, max_possible) ; o++){
                    if (seq2.compare(o-30, 30, end1) == 0){
                        if (closer(o, best_seed_hit)){
                            best_seed_hit = o;
                        }
                        if (closer(o, best_exact) && seq1.compare(seq1.size()-o, o, seq2, 0, o) == 0){
                            best_exact = o;
                        }
                    }
                }
            };
            int best_exact = -1;
            int best_seed_hit = -1;
            //first look in a window centered on the expected overlap: if an exact overlap is found there, the
            //closest exact overlap overall is necessarily in the window too
            long int window = std::max(200L, (long int) (0.15*expected_overlap));
            scan((long int) expected_overlap - window, (long int) expected_overlap + window, best_exact, best_seed_hit);
            if (best_exact == -1){
                best_seed_hit = -1;
                scan(30, max_possible, best_exact, best_seed_hit);
            }
            if (best_exact != -1){
                overlap = best_exact;
            }
            else if (best_seed_hit != -1){
                overlap = best_seed_hit;
            }
            else{
                overlap = min((long int) default_overlap, min(size1, size2));
            }
        }
        return "L\t" + name1 + "\t" + orientation1 + "\t" + name2 + "\t" + orientation2 + "\t" + std::to_string(overlap) + "M\n";
    };

    //compute the links by blocks of consecutive L lines, in parallel, and write them in order (the other lines, e.g. the
    //long S lines, are written directly and never buffered)
    const size_t BLOCK_SIZE = 100000;
    vector<string> block;
    vector<string> results;
    auto flush_block = [&](){
        if (block.empty()){
            return;
        }
        results.assign(block.size(), "");
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic, 64)
        for (size_t i = 0 ; i < block.size() ; i++){
            results[i] = compute_link(block[i], *sequence_readers[omp_get_thread_num()]);
        }
        for (const string& r : results){
            out << r;
        }
        block.clear();
    };
    while (std::getline(input, line)){
        if (line[0] == 'L'){
            block.push_back(line);
            if (block.size() == BLOCK_SIZE){
                flush_block();
            }
        }
        else{
            flush_block();
            out << line << "\n";
        }
    }
    flush_block();
}

/**
 * @brief inplace, put the contigs first and the links after
 * 
 * @param gfa 
 */
void sort_GFA(std::string gfa){
    ifstream input(gfa);
    ofstream output(gfa + ".sorted");
    string line;
    vector<string> contigs;
    vector<string> links;
    while (std::getline(input, line))
    {
        if (line[0] == 'S'){
            contigs.push_back(line);
        }
        else if (line[0] == 'L'){
            links.push_back(line);
        }
    }
    input.close();

    for (const auto& c: contigs){
        output << c << "\n";
    }
    for (const auto& l: links){
        output << l << "\n";
    }
    output.close();

    //move the sorted file to the original file
    if (std::rename((gfa + ".sorted").c_str(), gfa.c_str()) != 0){
        cerr << "ERROR: could not move " << gfa << ".sorted to " << gfa << "\n";
        exit(1);
    }
}


/**
 * @brief Equivalent to (forward ? seq : reverse_complement(seq)).substr(pos, len), without copying the whole sequence
 */
static string oriented_substr(const string& seq, bool forward, size_t pos, size_t len = string::npos){
    if (forward){
        return seq.substr(pos, len);
    }
    if (pos > seq.size()){
        throw std::out_of_range("oriented_substr");
    }
    size_t n = min(len, seq.size() - pos);
    string piece = seq.substr(seq.size() - pos - n, n);
    return reverse_complement(piece);
}

/**
 * @brief Graph with integer IDs for the contigs, shared by the functions that clean the graph and align reads on it
 * (instead of each of them building its own maps indexed by contig names)
 */
struct LinkedGraph {
    vector<string> names;
    unordered_flat_map<string, int> ids;
    vector<int> length; //0 for contigs only seen in L lines
    vector<string> sequences; //only if loaded
    vector<float> depth; //first DP or km tag of the S line
    vector<bool> has_depth;
    vector<bool> has_S_line;
    vector<std::array<vector<pair<int,char>>, 2>> links; //[0]: links of the left end, [1]: right end. Each link is (neighbor, end of the neighbor)
    vector<int> segments_in_order; //IDs in the order of the S lines

    int id_of(const string& name){
        auto it = ids.find(name);
        if (it != ids.end()){
            return it->second;
        }
        int id = names.size();
        ids[name] = id;
        names.push_back(name);
        length.push_back(0);
        sequences.push_back("");
        depth.push_back(0);
        has_depth.push_back(false);
        has_S_line.push_back(false);
        links.push_back({});
        return id;
    }
    int find(const string& name) const {
        auto it = ids.find(name);
        return it == ids.end() ? -1 : it->second;
    }
    int size() const {
        return names.size();
    }
};

static LinkedGraph load_linked_graph(const string& gfa_file, bool load_sequences){
    LinkedGraph graph;
    ifstream input(gfa_file);
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
            int id = graph.id_of(name);
            graph.length[id] = sequence.size();
            graph.has_S_line[id] = true;
            graph.segments_in_order.push_back(id);
            //first DP or km tag
            string tag;
            while (ss >> tag && tag.size() >= 2){
                if (tag.substr(0,2) == "DP" || tag.substr(0,2) == "km"){
                    graph.depth[id] = std::stof(tag.substr(5, tag.size()-5));
                    graph.has_depth[id] = true;
                    break;
                }
            }
            if (load_sequences){
                graph.sequences[id] = std::move(sequence);
            }
        }
        else if (line[0] == 'L'){
            string name1;
            string name2;
            string orientation1;
            string orientation2;
            string dont_care;
            std::stringstream ss(line);
            ss >> dont_care >> name1 >> orientation1 >> name2 >> orientation2;
            int id1 = graph.id_of(name1);
            int id2 = graph.id_of(name2);

            //end of each contig involved in the link
            char end1 = (orientation1 == "+" ? 1 : 0);
            char end2 = (orientation2 == "+" ? 0 : 1);
            pair<int,char> neighbor = {id2, end2};
            auto& links1 = graph.links[id1][end1];
            if (std::find(links1.begin(), links1.end(), neighbor) == links1.end()){
                links1.push_back(neighbor);
            }
            neighbor = {id1, end1};
            auto& links2 = graph.links[id2][end2];
            if (std::find(links2.begin(), links2.end(), neighbor) == links2.end()){
                links2.push_back(neighbor);
            }
        }
    }
    return graph;
}

//contigs must have a DP or km tag
static void check_depths(const LinkedGraph& graph){
    for (int id : graph.segments_in_order){
        if (!graph.has_depth[id]){
            cerr << "ERROR: no depth found for contig " << graph.names[id] << "\n";
            exit(1);
        }
    }
}

//links must be between contigs that have an S line
static void check_links(const LinkedGraph& graph){
    for (int id = 0 ; id < graph.size() ; id++){
        if (!graph.has_S_line[id]){
            cerr << "ERROR: contig not found in linked in pop_bubbles: " << graph.names[id] << "\n";
            exit(1);
        }
    }
}

//index every 5th kmer of the contigs: kmer -> (contig, position)
static unordered_flat_map<uint64_t, pair<int,int>> index_kmers_of_contigs(LinkedGraph& graph, int km){
    unordered_flat_map<uint64_t, pair<int,int>> kmers_to_contigs;
    for (int contig : graph.segments_in_order){
        uint64_t hash_foward = -1;
        size_t pos_end = 0;
        long pos_begin = -km;
        while (roll_f(hash_foward, km, graph.sequences[contig], pos_end, pos_begin, false)){
            if (pos_begin>=0  && pos_begin % 5 == 0){
                kmers_to_contigs[hash_foward] = make_pair(contig, pos_end-km);
            }
        }
    }
    return kmers_to_contigs;
}

static vector<vector<pair<int, bool>>> list_all_paths_from_contig(const LinkedGraph& graph, int start_contig, bool start_orientation, int max_length, int km, int target_contig, bool target_orientation);

struct Path{
    vector<int> contigs;
    vector<bool> orientations;
    int start_position_on_contig;
    int end_position_on_contig;
};

/**
 * @brief Create a gaf from unitig graph object and a set of reads. Also compute the coverage of the contigs
 * 
 * @param unitig_graph 
 * @param km 
 * @param reads_file 
 * @param output_file
 * @param hard_correct if true, only keep the reads that can be perfectly corrected (i.e. for which we can find a single path in the graph). If false, keep all reads but only correct the part of the read that can be unambiguously corrected
 * @param coverages
 */
void create_corrected_reads_from_unitig_graph(std::string unitig_graph, int km, std::string reads_file, std::string output_file, bool hard_correct, robin_hood::unordered_flat_map<std::string, float>& coverages, int num_threads){
    
    LinkedGraph graph = load_linked_graph(unitig_graph, true);
    unordered_flat_map<uint64_t, pair<int,int>> kmers_to_contigs = index_kmers_of_contigs(graph, km); //in what contig is the kmer and at what position (only unique kmer ofc, meant to work with unitig graph)
    vector<float> coverage_of_contigs (graph.size(), 0);

    // Prepare for parallel processing
    omp_lock_t coverage_lock;
    omp_init_lock(&coverage_lock);
    
    int nb_reads = 0;
    int nb_reads_single_path = 0;
    omp_lock_t count_lock;
    omp_init_lock(&count_lock);
    
    ofstream output(output_file);
    omp_lock_t output_lock;
    omp_init_lock(&output_lock);
    
    // Read file in chunks
    ifstream input2(reads_file);
    const int CHUNK_SIZE = 20000;
    bool done = false;
    
    while (!done) {
        // Read a chunk of reads (main thread only)
        vector<pair<string, string>> read_chunk; // (name, sequence)
        read_chunk.reserve(CHUNK_SIZE);
        
        string name_line, seq_line;
        for (int i = 0; i < CHUNK_SIZE && std::getline(input2, name_line); i++) {
            if (name_line.empty()) continue;
            
            if (name_line[0] == '@' || name_line[0] == '>') {
                string name = name_line.substr(1);
                if (name.substr(0, 5) == "false") continue;
                
                if (std::getline(input2, seq_line)) {
                    read_chunk.push_back({name, seq_line});
                }
            }
        }
        
        if (read_chunk.empty()) {
            done = true;
            break;
        }
        
        // Process chunk in parallel
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t read_idx = 0; read_idx < read_chunk.size(); read_idx++) {
            const string& name = read_chunk[read_idx].first;
            string line = read_chunk[read_idx].second; // Make a copy since roll() needs non-const reference
            
            string corrected_seq = "";
            int last_corrected_pos = 0; // Position up to which we've added to corrected_seq
            
            if (line.size() < km){
                continue;
            }

            uint64_t hash_foward = 0;
            uint64_t hash_reverse = 0;

            int pos_to_look_at = 0;
            size_t pos_end = 0;
            long pos_begin = -km;
            long pos_middle = -(km+1)/2;
            int previous_match_pos = -1;
            int previous_contig = -1;
            bool previous_orientation = true;
            int previous_match_pos_on_contig = -1;
            bool read_is_single_path = true;
            
            while(roll(hash_foward, hash_reverse, km, line, pos_end, pos_begin, pos_middle, false)){
                if (pos_begin == pos_to_look_at){

                    unsigned long kmer = hash_foward; 
                    bool found_match = false;
                    bool forward_orientation = true;

                    auto match = kmers_to_contigs.find(kmer);
                    if (match != kmers_to_contigs.end()){
                        found_match = true;
                        forward_orientation = true;
                    }
                    else {
                        match = kmers_to_contigs.find(hash_reverse);
                        if (match != kmers_to_contigs.end()){
                            found_match = true;
                            forward_orientation = false;
                            kmer = hash_reverse;
                        }
                    }

                    if (found_match){
                        int contig = match->second.first;
                        int pos_in_contig = match->second.second;
                        
                        if (!forward_orientation){
                            pos_in_contig = graph.length[contig] - pos_in_contig - km;
                        }

                        // Handle sequence before this match
                        if (previous_match_pos == -1){
                            // First match: add pos_begin bases of contig sequence
                            corrected_seq += oriented_substr(graph.sequences[contig], forward_orientation, std::max((long int) 0, pos_in_contig - pos_begin), min(pos_begin, (long int) pos_in_contig));
                            if (pos_in_contig - pos_begin < 0 && !hard_correct){
                                corrected_seq = line.substr(0, pos_begin - pos_in_contig) + corrected_seq;
                            }
                            last_corrected_pos = pos_begin;
                        }
                        else if (contig == previous_contig && forward_orientation == previous_orientation && pos_in_contig >= previous_match_pos_on_contig){ //do a special case because that is very frequent
                            // Add sequence between previous match and current match on the same contig
                            // Extract the portion between the two matches
                            corrected_seq += oriented_substr(graph.sequences[contig], forward_orientation, previous_match_pos_on_contig, pos_in_contig - previous_match_pos_on_contig);
                            last_corrected_pos = pos_begin;
                            // Update coverage for this portion
                            omp_set_lock(&coverage_lock);
                            coverage_of_contigs[contig] += (pos_in_contig - previous_match_pos_on_contig) / (double)graph.length[contig];
                            omp_unset_lock(&coverage_lock);
                        }
                        else{
                            // Try to bridge from previous match to current match
                            int distance = max((long int) 1,pos_begin - previous_match_pos - pos_in_contig); // is distance is negative, still look the neighboring contig

                            vector<vector<pair<int, bool>>> all_paths;
                            if (distance > 0){
                                all_paths = list_all_paths_from_contig(graph, contig, !forward_orientation, distance, km, previous_contig, !previous_orientation);
                            }
                            else{
                                all_paths = {{{contig, !forward_orientation}}};
                            }
                            
                            //the path must end on the previous contig, traversed in the orientation in which the read saw it
                            int valid_path_count = 0;
                            int valid_path_index = -1;
                            for (int path_idx = 0; path_idx < all_paths.size(); path_idx++) {
                                const auto& node = all_paths[path_idx][all_paths[path_idx].size()-1];
                                if (node.first == previous_contig && node.second == !previous_orientation && all_paths[path_idx].size() > 1) {
                                    valid_path_count++;
                                    valid_path_index = path_idx;
                                }
                            }
                            
                            if (valid_path_count == 1) {
                                // Found unique path - use bridging contigs
                                const auto& bridging_path = all_paths[valid_path_index];
                                string correct_seq = "";

                                // First append to the read the end of the contig that was mathched last (previous_contig)
                                const string& prev_seq = graph.sequences[previous_contig];
                                correct_seq += oriented_substr(prev_seq, previous_orientation, previous_match_pos_on_contig, prev_seq.size()-previous_match_pos_on_contig);
                                // Update coverage for the previous contig based on the portion used
                                omp_set_lock(&coverage_lock);
                                int length_used = min((long int) prev_seq.size()-previous_match_pos_on_contig, pos_begin-previous_match_pos+km-1);
                                coverage_of_contigs[previous_contig] += length_used / (double)graph.length[previous_contig];
                                omp_unset_lock(&coverage_lock);

                                for (int i = bridging_path.size() - 2; i >= 1; i--) {
                                    // Skip overlap
                                    correct_seq += oriented_substr(graph.sequences[bridging_path[i].first], !bridging_path[i].second, km-1);
                                    
                                    omp_set_lock(&coverage_lock);
                                    coverage_of_contigs[bridging_path[i].first] += 1;
                                    omp_unset_lock(&coverage_lock);
                                }

                                if (bridging_path.size() > 1){
                                    // Add the beginning of current contig to complete the bridging
                                    // Add up to pos_in_contig bases from the current contig
                                    correct_seq += oriented_substr(graph.sequences[contig], forward_orientation, km-1, pos_in_contig);
                                    omp_set_lock(&coverage_lock);
                                    int length_used = pos_in_contig;
                                    coverage_of_contigs[contig] += length_used / (double)graph.length[contig];
                                    omp_unset_lock(&coverage_lock);
                                }
                                correct_seq = correct_seq.substr(0, correct_seq.size()-km+1);
                                corrected_seq += correct_seq;

                                last_corrected_pos = pos_begin;
                            }
                            else {        
                                if (!hard_correct){
                                    corrected_seq += line.substr(last_corrected_pos, pos_begin - last_corrected_pos);
                                }
                                else {
                                    //output the corrected seq we have and reset the corrected seq for the next path
                                    if (corrected_seq.size() > 0){
                                        stringstream thread_output;
                                        thread_output << ">" << name << "\n";
                                        thread_output << corrected_seq << "\n";
                                        // Write to output file with lock
                                        omp_set_lock(&output_lock);
                                        output << thread_output.str();
                                        omp_unset_lock(&output_lock);
                                    }
                                    corrected_seq = "";
                                }
                                last_corrected_pos = pos_begin;
                                read_is_single_path = false;
                            }
                        }

                        int length_of_contig_left = graph.length[contig] - pos_in_contig - km;
                        if (length_of_contig_left > 10){
                            pos_to_look_at += min((int)(length_of_contig_left*0.8), (int)(line.size() - pos_begin - km - 5));
                        }

                        previous_match_pos = pos_begin;
                        previous_contig = contig;
                        previous_orientation = forward_orientation;
                        previous_match_pos_on_contig = pos_in_contig;
                    }
                    pos_to_look_at++;
                }
            }

            // Add any remaining sequence at the end
            if (last_corrected_pos < line.size()){
                // Add the sequence of the contig after the last match if it is long enough. If not, add the remaining read sequence (only if not in hard_correct mode), or add the rest of the contig (if in hard_correct mode)
                if (previous_match_pos != -1){
                    int remaining_contig_length = graph.length[previous_contig] - previous_match_pos_on_contig;
                    int remaining_read_length = line.size() - last_corrected_pos;
                    corrected_seq += oriented_substr(graph.sequences[previous_contig], previous_orientation, previous_match_pos_on_contig, min(remaining_contig_length, remaining_read_length));
                    if (remaining_contig_length < remaining_read_length && !hard_correct){
                        corrected_seq += line.substr(last_corrected_pos + remaining_contig_length);
                    }
                }
                else{
                    read_is_single_path = false;
                    if (!hard_correct){
                        corrected_seq += line.substr(last_corrected_pos);
                    }
                } 
            }

            // Update counters with lock
            omp_set_lock(&count_lock);
            nb_reads++;
            if (read_is_single_path){
                nb_reads_single_path++;
            }
            omp_unset_lock(&count_lock);

            // Build output string for this read
            stringstream thread_output;
            if (corrected_seq.size() > 0){
                // FASTA format output (corrected read)
                thread_output << ">" << name << "\n";
                thread_output << corrected_seq << "\n";
            }
            // Write to output file with lock
            if (thread_output.str().size() > 0) {
                omp_set_lock(&output_lock);
                output << thread_output.str();
                omp_unset_lock(&output_lock);
            }
        } // end parallel for
    } // end while chunks

    input2.close();
    output.close();
    
    omp_destroy_lock(&coverage_lock);
    omp_destroy_lock(&count_lock);
    omp_destroy_lock(&output_lock);

    for (int contig = 0 ; contig < graph.size() ; contig++){
        coverages[graph.names[contig]] = coverage_of_contigs[contig];
    }
    add_coverages_to_graph(unitig_graph, coverages);

    cout << "    -> Number of cleanly corrected/aligned reads: " << nb_reads_single_path << " out of " << nb_reads << endl;
}

void create_gaf_from_unitig_graph(std::string unitig_graph, int km, std::string reads_file, std::string output_file, robin_hood::unordered_flat_map<std::string, float>& coverages, int num_threads){
    
    LinkedGraph graph = load_linked_graph(unitig_graph, true);
    unordered_flat_map<uint64_t, pair<int,int>> kmers_to_contigs = index_kmers_of_contigs(graph, km);
    vector<float> coverage_of_contigs (graph.size(), 0);

    omp_lock_t coverage_lock;
    omp_init_lock(&coverage_lock);
    
    int nb_reads = 0;
    int nb_reads_single_path = 0;
    omp_lock_t count_lock;
    omp_init_lock(&count_lock);
    
    ofstream output(output_file);
    omp_lock_t output_lock;
    omp_init_lock(&output_lock);
    
    ifstream input2(reads_file);
    const int CHUNK_SIZE = 20000;
    bool done = false;
    
    while (!done) {
        vector<pair<string, string>> read_chunk;
        read_chunk.reserve(CHUNK_SIZE);
        
        string name_line, seq_line;
        for (int i = 0; i < CHUNK_SIZE && std::getline(input2, name_line); i++) {
            if (name_line.empty()) continue;
            
            if (name_line[0] == '@' || name_line[0] == '>') {
                string name = name_line.substr(1);
                if (name.substr(0, 5) == "false") continue;
                
                if (std::getline(input2, seq_line)) {
                    read_chunk.push_back({name, seq_line});
                }
            }
        }
        
        if (read_chunk.empty()) {
            done = true;
            break;
        }
        
        #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
        for (size_t read_idx = 0; read_idx < read_chunk.size(); read_idx++) {
            const string& name = read_chunk[read_idx].first;
            string line = read_chunk[read_idx].second;
            
            vector<Path> paths;
            Path current_path;
            current_path.start_position_on_contig = -1;
            current_path.end_position_on_contig = -1;
            
            if (line.size() < km){
                continue;
            }

            uint64_t hash_foward = 0;
            uint64_t hash_reverse = 0;

            int pos_to_look_at = 0;
            size_t pos_end = 0;
            long pos_begin = -km;
            long pos_middle = -(km+1)/2;
            int previous_match = 0;
            int previous_contig = -1;
            bool previous_forward = true; //orientation in which the read saw previous_contig
            
            while(roll(hash_foward, hash_reverse, km, line, pos_end, pos_begin, pos_middle, false)){
                if (pos_begin == pos_to_look_at){

                    unsigned long kmer = hash_foward; 

                    auto match_fw = kmers_to_contigs.find(kmer);
                    auto match_rv = match_fw == kmers_to_contigs.end() ? kmers_to_contigs.find(hash_reverse) : kmers_to_contigs.end();
                    if (match_fw != kmers_to_contigs.end()){
                        int contig = match_fw->second.first;
                        int pos_in_contig = match_fw->second.second;

                        if (current_path.contigs.size() == 0){
                            current_path.start_position_on_contig = pos_in_contig;
                        }
                        
                        //same contig, moving forward on it: the read is still on it. If the position moved backward, the read went around a loop (e.g. a circular contig)
                        if (current_path.contigs.size() > 0 && current_path.contigs[current_path.contigs.size()-1] == contig && current_path.orientations[current_path.orientations.size()-1] == true
                            && pos_in_contig + km >= current_path.end_position_on_contig){
                            current_path.end_position_on_contig = pos_in_contig + km;
                        }
                        else{
                            if (previous_contig != -1){
                                int distance = pos_to_look_at - previous_match - pos_in_contig;
                                vector<vector<pair<int, bool>>> all_paths = list_all_paths_from_contig(graph, contig, false, distance, km, previous_contig, !previous_forward);
                                
                                int valid_path_count = 0;
                                int valid_path_index = -1;
                                for (int path_idx = 0; path_idx < all_paths.size(); path_idx++) {
                                    const auto& node = all_paths[path_idx][all_paths[path_idx].size()-1];
                                    if (node.first == previous_contig && node.second == !previous_forward) {
                                        valid_path_count++;
                                        valid_path_index = path_idx;
                                    }
                                }
                                
                                if (valid_path_count == 1) {
                                    const auto& path = all_paths[valid_path_index];
                                    for (int i = path.size() - 2; i >= 0; i--) {
                                        if (path[i].first != contig) {
                                            current_path.contigs.push_back(path[i].first);
                                            current_path.orientations.push_back(!path[i].second);
                                            omp_set_lock(&coverage_lock);
                                            coverage_of_contigs[path[i].first] += min(1.0, (line.size() - pos_begin) / (double)graph.length[path[i].first]);
                                            omp_unset_lock(&coverage_lock);
                                        }
                                    }
                                }
                                
                                if (valid_path_count != 1 && current_path.contigs.size() > 0) {
                                    paths.push_back(current_path);
                                    current_path.contigs.clear();
                                    current_path.orientations.clear();
                                    current_path.start_position_on_contig = pos_in_contig;
                                    current_path.end_position_on_contig = pos_in_contig + km;
                                }
                            }

                            current_path.contigs.push_back(contig);
                            current_path.orientations.push_back(true);
                            if (current_path.start_position_on_contig == -1) {
                                current_path.start_position_on_contig = pos_in_contig;
                            }
                            current_path.end_position_on_contig = pos_in_contig + km;
                            
                            omp_set_lock(&coverage_lock);
                            coverage_of_contigs[contig] += std::min(1.0, (line.size()- pos_begin) / (double)graph.length[contig]);
                            omp_unset_lock(&coverage_lock);
                        }
                        
                        int length_of_contig_left = graph.length[contig] - pos_in_contig - km;
                        if (length_of_contig_left > 10){
                            pos_to_look_at += min((int) (length_of_contig_left*0.8) , (int)(line.size() - pos_begin - km - 5));
                        }
                        previous_match = pos_begin;
                        previous_contig = contig;
                        previous_forward = true;
                    }
                    else if (match_rv != kmers_to_contigs.end()){
                        int contig = match_rv->second.first;
                        int pos_in_contig = match_rv->second.second;

                        if (current_path.contigs.size() == 0){
                            current_path.start_position_on_contig = graph.length[contig] - pos_in_contig - km;
                        }
                        
                        if (current_path.contigs.size() > 0 && current_path.contigs[current_path.contigs.size()-1] == contig && current_path.orientations[current_path.orientations.size()-1] == false
                            && graph.length[contig] - pos_in_contig >= current_path.end_position_on_contig){
                            current_path.end_position_on_contig = graph.length[contig] - pos_in_contig;
                        }
                        else{
                            if (previous_contig != -1){
                                int distance = pos_to_look_at - previous_match - (graph.length[contig] - pos_in_contig - km);
                                vector<vector<pair<int, bool>>> all_paths = list_all_paths_from_contig(graph, contig, true, distance, km, previous_contig, !previous_forward);

                                int valid_path_count = 0;
                                int valid_path_index = -1;
                                for (int path_idx = 0; path_idx < all_paths.size(); path_idx++) {
                                    const auto& node = all_paths[path_idx][all_paths[path_idx].size()-1];
                                    if (node.first == previous_contig && node.second == !previous_forward) {
                                        valid_path_count++;
                                        valid_path_index = path_idx;
                                    }
                                }
                                
                                if (valid_path_count == 1) {
                                    const auto& path = all_paths[valid_path_index];
                                    for (int i = path.size() - 2; i >= 0; i--) {
                                        if (path[i].first != contig) {
                                            current_path.contigs.push_back(path[i].first);
                                            current_path.orientations.push_back(!path[i].second);
                                            omp_set_lock(&coverage_lock);
                                            coverage_of_contigs[path[i].first] += min(1.0, (line.size() - pos_begin) / (double)graph.length[path[i].first]);
                                            omp_unset_lock(&coverage_lock);
                                        }
                                    }
                                }
                                
                                if (valid_path_count != 1 && current_path.contigs.size() > 0) {
                                    paths.push_back(current_path);
                                    current_path.contigs.clear();
                                    current_path.orientations.clear();
                                    current_path.start_position_on_contig = graph.length[contig] - pos_in_contig - km;
                                    current_path.end_position_on_contig = graph.length[contig] - pos_in_contig;
                                }
                            }

                            current_path.contigs.push_back(contig);
                            current_path.orientations.push_back(false);
                            if (current_path.start_position_on_contig == -1) {
                                current_path.start_position_on_contig = graph.length[contig] - pos_in_contig - km;
                            }
                            current_path.end_position_on_contig = graph.length[contig] - pos_in_contig;
                            
                            omp_set_lock(&coverage_lock);
                            coverage_of_contigs[contig] += min(1.0, (line.size()- pos_begin) / (double)graph.length[contig]);
                            omp_unset_lock(&coverage_lock);
                        }
                        
                        int length_of_contig_left = pos_in_contig;
                        if (length_of_contig_left > 10){
                            pos_to_look_at += min((int) (length_of_contig_left*0.8), (int)(line.size() - pos_begin - km - 5));
                        }
                        previous_match = pos_begin;
                        previous_contig = contig;
                        previous_forward = false;
                    }
                    pos_to_look_at++;
                }
            }

            if (current_path.contigs.size() > 0){
                paths.push_back(current_path);
            }

            omp_set_lock(&count_lock);
            nb_reads++;
            if (paths.size() == 1){
                nb_reads_single_path++;
            }
            omp_unset_lock(&count_lock);

            stringstream thread_output;
            int idx_of_path = 0;
            for (const auto& p : paths) {
                if (p.contigs.size() > 0){
                    thread_output << name << "_" << idx_of_path << "\t" << line.size() << "\t0\t" << line.size() << "\t+\t";
                    
                    for (int i = 0; i < p.contigs.size(); i++){
                        thread_output << (p.orientations[i] ? ">" : "<") << graph.names[p.contigs[i]];
                    }
                    thread_output << "\t";
                    
                    int path_length = 0;
                    for (const auto& contig : p.contigs){
                        path_length += graph.length[contig];
                    }
                    thread_output << path_length << "\t0\t" << path_length << "\t" << line.size() << "\t" << line.size() << "\t255\n";
                    
                    idx_of_path += 1;
                }
            }
            
            if (thread_output.str().size() > 0) {
                omp_set_lock(&output_lock);
                output << thread_output.str();
                omp_unset_lock(&output_lock);
            }
        }
    }

    input2.close();
    output.close();
    
    omp_destroy_lock(&coverage_lock);
    omp_destroy_lock(&count_lock);
    omp_destroy_lock(&output_lock);

    for (int contig = 0 ; contig < graph.size() ; contig++){
        coverages[graph.names[contig]] = coverage_of_contigs[contig];
    }
    add_coverages_to_graph(unitig_graph, coverages);

    cout << "    -> Number of cleanly corrected/aligned reads: " << nb_reads_single_path << " out of " << nb_reads << endl;
}


/**
 * @brief Given a starting position on a contig and an orientation, follow the graph and list all possible contigs and paths
 * 
 * @param graph
 * @param start_contig starting contig
 * @param start_orientation true if we start from the right end ('+'), false if from the left end ('-')
 * @param max_length maximum total length of contigs to explore in the paths
 * @param target_contig if not -1, a path stops as soon as it reaches this contig in target_orientation
 * @return vector<vector<pair<int, bool>>> list of paths, each path is a list of (contig, orientation)
 */
static vector<vector<pair<int, bool>>> list_all_paths_from_contig(const LinkedGraph& graph, int start_contig, bool start_orientation, int max_length, int km, int target_contig, bool target_orientation){
    
    vector<vector<pair<int, bool>>> all_results;
    all_results.reserve(50); // Reserve space for expected maximum paths
    
    // Recursive exploration function using push/pop
    std::function<void(int, char, int, vector<pair<int, bool>>&)> explore;
    explore = [&](int current_contig, char current_end, int length_left, vector<pair<int, bool>>& current_path) {
        
        if (length_left <= 0){
            all_results.push_back(current_path);
            return;
        }
        
        // Early termination if too many paths
        if (all_results.size() >= 50){
            return;
        }
        
        // Get neighbors from the appropriate end
        const vector<pair<int, char>>& neighbors = graph.links[current_contig][current_end == 1 ? 1 : 0];
        
        if (neighbors.size() == 0){
            all_results.push_back(current_path);
            return;
        }
        
        for (const auto& neighbor : neighbors){
            // neighbor.second tells us which end of the neighbor we arrive at
            // if we arrive at end 0, we traverse it in forward orientation (true)
            // if we arrive at end 1, we traverse it in reverse orientation (false)
            bool neighbor_orientation = (neighbor.second == 0);
            
            // Push to path
            current_path.push_back({neighbor.first, neighbor_orientation});
            
            if (target_contig != -1 && neighbor.first == target_contig && neighbor_orientation == target_orientation){
                //reached the target: do not look further, even if the length budget is not exactly exhausted (indels in the reads)
                all_results.push_back(current_path);
            }
            else{
                // Continue exploration from the opposite end of the neighbor
                char next_end = 1 - neighbor.second;
                //at least 1, so that the exploration ends even in cycles of contigs shorter than k (or of unknown length)
                int neighbor_length = max(1, graph.length[neighbor.first] - km + 1);
                
                explore(neighbor.first, next_end, length_left - neighbor_length, current_path);
            }
            
            // Pop from path (backtrack)
            current_path.pop_back();
            
            // Early termination check
            if (all_results.size() >= 50){
                return;
            }
        }
    };
    
    // Start exploration
    vector<pair<int, bool>> initial_path;
    initial_path.reserve(20); // Reserve space for typical path length
    initial_path.push_back({start_contig, start_orientation});
    
    char start_end = start_orientation ? 1 : 0; // if orientation is true ('+'), we start from end 1 (right)
    
    explore(start_contig, start_end, max_length, initial_path);
    
    // Return empty if we hit the limit
    if (all_results.size() >= 50){
        return {};
    }
    
    return all_results;
}

void gfa_to_fasta(string gfa, string fasta){   
        ifstream input(gfa);
        ofstream out(fasta);
    
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
                
                out << ">" << name << "\n";
                out << sequence << "\n";
            }
        }
}

/**
 * @brief Given a gfa file and coverage of contigs, append the dp:f: tag to the S lines of the gfa file, suppressing other dp tags or kc tags
 * 
 * @param gfa 
 * @param coverages 
 */
void add_coverages_to_graph(std::string gfa, robin_hood::unordered_map<std::string, float>& coverages){

    ifstream input(gfa);
    ofstream out(gfa + ".tmp");
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
            float coverage = 0;
            if (coverages.find(name) != coverages.end()){
                coverage = coverages[name];
            }
            out << "S\t" << name << "\t" << sequence << "\tDP:f:" << coverage << "\tLN:i:" << sequence.size() << "\n";
        }
        else{
            out << line << "\n";
        }
    }
    out.close();

    //move the file with coverages to the original file
    if (std::rename((gfa + ".tmp").c_str(), gfa.c_str()) != 0){
        cerr << "ERROR: could not move " << gfa << ".tmp to " << gfa << "\n";
        exit(1);
    }
}

//flags of the memo of pop_and_shave_graph: a path that is not overcovered was found left / right of the contig
static const uint8_t NOT_OVERCOVERED_LEFT = 1;
static const uint8_t NOT_OVERCOVERED_RIGHT = 2;

/**
 * @brief Recursive exploration of one side of a contig, within length_left bases, looking for a path that does not go
 * through a contig with a coverage above big_coverage
 * 
 * @param contig current contig
 * @param endOfContig endOfContig we arrive from
 * @param links links of the graph (see LinkedGraph)
 * @param coverage coverage of each contig
 * @param length_of_contigs length of each contig
 * @param k 
 * @param length_left 
 * @param original_coverage coverage of the contig whose neighborhood is explored
 * @param big_coverage coverage above which a contig is overcovered
 * @param not_overcovered flags (NOT_OVERCOVERED_LEFT/RIGHT) of the contigs already known not to be overcovered left or right
 * @param memo results of the states already explored during this exploration
 * @return int 2: not overcovered path found, 1: overcovered, 0: nothing overcovered but dead end
 */
static int explore_neighborhood_uncached(int contig, int endOfContig, const vector<std::array<vector<pair<int,char>>, 2>>& links, const vector<float>& coverage, const vector<int>& length_of_contigs, int k, int length_left, double original_coverage, double big_coverage, vector<std::atomic<uint8_t>>& not_overcovered, unordered_flat_map<uint64_t, int8_t>& memo);

/**
 * @brief Memoized exploration of the neighborhood: within one exploration (fixed original_coverage and big_coverage),
 * the result only depends on (contig, end, length left), so each state is explored once instead of once per path
 * leading to it (which is exponential in tangles of short contigs)
 */
static int explore_neighborhood(int contig, int endOfContig, const vector<std::array<vector<pair<int,char>>, 2>>& links, const vector<float>& coverage, const vector<int>& length_of_contigs, int k, int length_left, double original_coverage, double big_coverage, vector<std::atomic<uint8_t>>& not_overcovered, unordered_flat_map<uint64_t, int8_t>& memo){
    uint64_t key = ((uint64_t) contig << 33) | ((uint64_t) (endOfContig & 1) << 32) | (uint32_t) length_left;
    auto it = memo.find(key);
    if (it != memo.end()){
        return it->second;
    }
    int result = explore_neighborhood_uncached(contig, endOfContig, links, coverage, length_of_contigs, k, length_left, original_coverage, big_coverage, not_overcovered, memo);
    memo[key] = result;
    return result;
}

static int explore_neighborhood_uncached(int contig, int endOfContig, const vector<std::array<vector<pair<int,char>>, 2>>& links, const vector<float>& coverage, const vector<int>& length_of_contigs, int k, int length_left, double original_coverage, double big_coverage, vector<std::atomic<uint8_t>>& not_overcovered, unordered_flat_map<uint64_t, int8_t>& memo){
    
    if (length_left <= 0){

        if (coverage[contig] > big_coverage){
            return 1;
        }
        else if (endOfContig == 1 && links[contig][0].size() == 0 || endOfContig == 0 && links[contig][1].size() == 0){
            return 0;
        }
        else{
            return 2;
        }
    }

    bool overcovered_in_neighborhood = false;
    if (endOfContig == 1){ //arrived by the right, go through to the left
        if (links[contig][0].size() == 0){
            return 0;
        }
        for (const auto& l: links[contig][0]){
            if (coverage[l.first] > big_coverage){
                overcovered_in_neighborhood = true;
            }
            else if (l.first != contig) { //the condition is so we don't end up in a loop

                //check if we're encountering an already known not overcovered contig
                uint8_t known = not_overcovered[l.first];
                if ((known & NOT_OVERCOVERED_LEFT) && l.second == 1 && coverage[l.first] < 2*original_coverage && coverage[l.first]*2 > original_coverage){
                    return 2;
                }
                else if ((known & NOT_OVERCOVERED_RIGHT) && l.second == 0 && coverage[l.first] < 2*original_coverage && coverage[l.first]*2 > original_coverage){
                    return 2;
                }

                int res = explore_neighborhood(l.first, l.second, links, coverage, length_of_contigs, k, length_left - max(1, length_of_contigs[l.first] - k + 1), original_coverage, big_coverage, not_overcovered, memo);
                if (res == 2){
                    return 2;
                }
                if (res == 1){
                    overcovered_in_neighborhood = true;
                }
            }
        }
        if (overcovered_in_neighborhood){
            return 1;
        }
        else{
            return 0;
        }
    }
    else if (endOfContig == 0){ //arrived by the left, go through to the right
        if (links[contig][1].size() == 0){
            return 0;
        }
        for (const auto& l: links[contig][1]){
            if (coverage[l.first] > big_coverage){
                overcovered_in_neighborhood = true;
            }
            else if (l.first != contig) { //the condition is so we don't end up in a loop

                uint8_t known = not_overcovered[l.first];
                if ((known & NOT_OVERCOVERED_LEFT) && l.second == 1 && coverage[l.first] < 2*original_coverage && coverage[l.first]*2 > original_coverage){
                    return 2;
                }
                else if ((known & NOT_OVERCOVERED_RIGHT) && l.second == 0 && coverage[l.first] < 2*original_coverage && coverage[l.first]*2 > original_coverage){
                    return 2;
                }

                int res = explore_neighborhood(l.first, l.second, links, coverage, length_of_contigs, k, length_left - max(1, length_of_contigs[l.first] - k + 1), original_coverage, big_coverage, not_overcovered, memo);
                if (res == 2){
                    return 2;
                }
                if (res == 1){
                    overcovered_in_neighborhood = true;
                }
            }
        }
        if (overcovered_in_neighborhood){
            return 1;
        }
        else{
            return 0;
        }
    }
    else{
        cerr << "ERROR: endOfContig should be 0 or 1 in explore_neighborhood\n";
        exit(1);
    }
    
}

/**
 * @brief Function that takes as input a graph, trim the tips less frequent than abundance_min, remove the branch of bubbles less abundant than abundance_min
 * 
 * @param gfa_in 
 * @param abundance_min //contigs above this coverage are solid, if -1, then coverage alone cannot save a contig
 * @param min_length  //contigs above this length are solid
 * @param contiguity //to collapse more bubbles
 * @param k
 * @param gfa_out 
 * @param extra_coverage //to retreat to the coverage because it comes from extra contigs added to the reads from previous assembly rounds
 */
void pop_and_shave_graph(string gfa_in, int abundance_min, int min_length, int k, string gfa_out, int extra_coverage, int num_threads, bool single_genome){

    if (min_length == -1){ //then no length is sufficient to keep a contig, it has to be done on the coverage
        //set min_length to the max int
        min_length = std::numeric_limits<int>::max();
    }

    LinkedGraph graph = load_linked_graph(gfa_in, false);
    check_depths(graph);
    const vector<std::array<vector<pair<int,char>>, 2>>& links = graph.links;
    const vector<int>& length_of_contigs = graph.length;
    const vector<int>& list_of_contigs = graph.segments_in_order;
    vector<float> coverage (graph.size(), 0);
    for (int contig : list_of_contigs){
        double depth = graph.depth[contig];
        coverage[contig] = std::max(1.0, depth - extra_coverage);
    }

    //kept contigs are the one with a coverage above abundance_min and their necessary neighbors for the contiguity
    vector<std::atomic<bool>> to_keep (graph.size());
    for (auto& k : to_keep){
        k = false;
    }

    //iterative cleaning of the graph

    vector<std::atomic<uint8_t>> not_overcovered (graph.size()); //iteratively mark the contigs that are not overcovered left or right
    for (auto& n : not_overcovered){
        n = 0;
    }

    //decide which contig we really want to keep
    #pragma omp parallel for num_threads(num_threads)
    for (int c = 0 ; c < list_of_contigs.size() ; c++){

        int contig = list_of_contigs[c];

        //check if this is a badly covered bubble
        bool bubble = false;
        if (links[contig][0].size() == 1 && links[contig][1].size() == 1){

            int neighbor_left = links[contig][0][0].first;
            char end_of_neighbor_left = links[contig][0][0].second;
            int neighbor_right = links[contig][1][0].first;
            char end_of_neighbor_right = links[contig][1][0].second;

            int other_neighbor_of_contig_left = -1;
            if (end_of_neighbor_left == 0 && links[neighbor_left][0].size() == 2){
                for (const auto& l: links[neighbor_left][0]){
                    if (l.first != contig){
                        other_neighbor_of_contig_left = l.first;
                    }
                }
            }
            else if (end_of_neighbor_left == 1 && links[neighbor_left][1].size() == 2){
                for (const auto& l: links[neighbor_left][1]){
                    if (l.first != contig){
                        other_neighbor_of_contig_left = l.first;
                    }
                }
            }

            int other_neighbor_of_contig_right = -1;
            if (end_of_neighbor_right == 0 && links[neighbor_right][0].size() == 2){
                for (const auto& l: links[neighbor_right][0]){
                    if (l.first != contig){
                        other_neighbor_of_contig_right = l.first;
                    }
                }
            }
            else if (end_of_neighbor_right == 1 && links[neighbor_right][1].size() == 2){
                for (const auto& l: links[neighbor_right][1]){
                    if (l.first != contig){
                        other_neighbor_of_contig_right = l.first;
                    }
                }
            }

            if (other_neighbor_of_contig_left == other_neighbor_of_contig_right && other_neighbor_of_contig_left != -1 && 5*coverage[contig] < coverage[other_neighbor_of_contig_left]){
                bubble = true;
            }
            if (single_genome && other_neighbor_of_contig_left == other_neighbor_of_contig_right && other_neighbor_of_contig_left != -1 && 2*coverage[contig] < coverage[other_neighbor_of_contig_left]){
                bubble = true;
            }
    
        }

        if (bubble && ((coverage[contig] < abundance_min || abundance_min == -1) && length_of_contigs[contig] < min_length)) //pop it if 1) this is a bubble and 2) coverage less than 5x and the contig is shorter than min_length
        {
            //do nothing, and most importantly, do not add the contig to the to_keep set
        }
        else if ((abundance_min != -1 && coverage[contig] > abundance_min) || length_of_contigs[contig] > min_length){ 
            
            int size_of_neighborhood = 7*k;
            //coverage above which a neighbor makes this contig look like an error. The better covered the contig, the more
            //overwhelming the neighbor must be: well-covered contigs squeezed between copies of a high-copy repeat are real
            double big_coverage = 200*coverage[contig];
            if (coverage[contig] < 5){ //not very solid, don't make such a fuss about deleting it
                big_coverage = 3*coverage[contig];
            }
            else if (coverage[contig] <= 10){
                big_coverage = 20*coverage[contig];
            }
            else if (coverage[contig] <= 20){
                big_coverage = 50*coverage[contig];
            }
            // cout << "launching..\n";
            int overcovered_right = 2;
            std::atomic<uint8_t>& not_overcovered_contig = not_overcovered[contig];
            if (!(not_overcovered_contig & NOT_OVERCOVERED_RIGHT)){
                unordered_flat_map<uint64_t, int8_t> memo;
                overcovered_right = explore_neighborhood(contig, 0, links, coverage, length_of_contigs, k, size_of_neighborhood, coverage[contig], big_coverage, not_overcovered, memo);
            }
            if (overcovered_right == 2){
                not_overcovered_contig |= NOT_OVERCOVERED_RIGHT;
            }
            int overcovered_left = 2;
            if (!(not_overcovered_contig & NOT_OVERCOVERED_LEFT)){
                unordered_flat_map<uint64_t, int8_t> memo;
                overcovered_left = explore_neighborhood(contig, 1, links, coverage, length_of_contigs, k, size_of_neighborhood, coverage[contig], big_coverage, not_overcovered, memo);
            }
            if (overcovered_left == 2){
                not_overcovered_contig |= NOT_OVERCOVERED_LEFT;
            }

            //now decide if the contig is to be kept
            if (overcovered_left == 1 && overcovered_right == 1){ //overcovered on both sides, pop it

            }
            else if (overcovered_left == 0 && overcovered_right == 0 && length_of_contigs[contig] > min_length){ //if the contig has two dead ends, keep it under conditiosn that it is long enough
                to_keep[contig] = true;
            }
            else if (overcovered_left == 2 || overcovered_right == 2){ //if the contig is not overcovered on one side, keep it (and it also passed the abundance_min or the min_length threshold)
                to_keep[contig] = true;
            }
            else if (overcovered_left == 1 && overcovered_right == 0 || overcovered_left == 0 && overcovered_right == 1){ //this means that this is a tip
                //do nothing, and most importantly, do not add the contig to the to_keep set
            }

            //if single genome mode is on, delete badly covered dead ends
            if (single_genome) {
                if (links[contig][0].size() == 0 && links[contig][1].size() > 0  ){
                    float max_neighbor_coverage = 0;
                    for (auto neighbor : links[contig][1]) {
                        max_neighbor_coverage = std::max(max_neighbor_coverage, coverage[neighbor.first]);
                    }
                    if (coverage[contig] * 2 < max_neighbor_coverage && length_of_contigs[contig] < 2*(long long)min_length) {
                        to_keep[contig] = false;
                    }
                }
                if (links[contig][0].size() > 0 && links[contig][1].size() == 0  ){
                    float max_neighbor_coverage = 0;
                    for (auto neighbor : links[contig][0]) {
                        max_neighbor_coverage = std::max(max_neighbor_coverage, coverage[neighbor.first]);
                    }
                    if (coverage[contig] * 2 < max_neighbor_coverage && length_of_contigs[contig] < 2*(long long)min_length) {
                        to_keep[contig] = false;
                    }
                }
                
            }
        }
    }

    //now add contigs that are necessary for the contiguity
    int number_of_edits = 1;
    while (number_of_edits > 0){
        number_of_edits = 0;

        #pragma omp parallel for num_threads(num_threads) reduction(+:number_of_edits)
        for (int c = 0 ; c < list_of_contigs.size() ; c++){

            int contig = list_of_contigs[c];

            if (to_keep[contig]){
                //make sure the contig has at least one neighbor left and right (if not, take the one with the highest coverage) 
                //this is the only way to keep contigs that are below abundance_min
                float best_coverage = 0;
                int best_contig = -1;
                bool at_least_one_neighbor = false;
                for (const auto& l: links[contig][1]){
                    if (coverage[l.first] > best_coverage){
                        best_coverage = coverage[l.first];
                        best_contig = l.first;
                    }
                    if (to_keep[l.first]){
                        at_least_one_neighbor = true;
                    }
                }
                if (best_contig != -1 && !at_least_one_neighbor){
                    if (!to_keep[best_contig].exchange(true)){
                        number_of_edits++;
                    }
                }

                best_coverage = 0;
                best_contig = -1;
                at_least_one_neighbor = false;
                for (const auto& l: links[contig][0]){
                    if (coverage[l.first] > best_coverage){
                        best_coverage = coverage[l.first];
                        best_contig = l.first;
                    }
                    if (to_keep[l.first]){
                        at_least_one_neighbor = true;
                    }
                }
                if (best_contig != -1 && !at_least_one_neighbor){
                    if (!to_keep[best_contig].exchange(true)){
                        number_of_edits++;
                    }
                }
            }
        }
    }

    //now write the gfa file without the contigs to remove
    ifstream input(gfa_in);
    string line;
    ofstream out(gfa_out);
    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            int id = graph.find(name);
            if (to_keep[id]){
                out << "S\t" << name << "\t" << sequence << "\tLN:i:" << sequence.size() << "\tkm:f:" << coverage[id] << "\n";
            }
        }
        else if (line[0] == 'L'){
            string name1;
            string name2;
            string orientation1;
            string orientation2;
            string dont_care;
            std::stringstream ss(line);
            int i = 0;
            ss >> dont_care >> name1 >> orientation1 >> name2 >> orientation2;

            if (to_keep[graph.find(name1)] && to_keep[graph.find(name2)]){
                out << line << "\n";
            }
        }
    }
    out.close();
}

/**
 * @brief At the ends of contigs that have several neighbors, cut the links to the neighbors that look like errors
 * (much less covered than the other neighbors, or dead ends that look like error tips)
 * 
 * @param gfa_in 
 * @param gfa_out 
 * @param k k of the graph (lengths are in compressed bases)
 */
void cut_links_for_contiguity(std::string gfa_in, std::string gfa_out, int k){

    LinkedGraph graph = load_linked_graph(gfa_in, false);
    check_depths(graph);
    check_links(graph);
    const vector<std::array<vector<pair<int,char>>, 2>>& links = graph.links;
    const vector<float>& coverage = graph.depth;
    const vector<int>& length_of_contigs = graph.length;

    //now cut the links that are not good for contiguity

    std::set<pair<pair<int,char>,pair<int,char>>> links_to_delete;
    for (int contig_name = 0 ; contig_name < graph.size() ; contig_name++){
        for (char end = 0 ; end < 2 ; end++){
            const vector<pair<int, char>>& neighbors = (end == 1) ? links[contig_name][1] : links[contig_name][0];
            if (neighbors.size() > 1){
                float max_coverage = 0;
                bool all_neighbors_are_dead_ends = true;
                for (auto& neighbor : neighbors){
                    if (coverage[neighbor.first] > max_coverage){
                        max_coverage = coverage[neighbor.first];
                    }
                    if (links[neighbor.first][0].size() > 0 && links[neighbor.first][1].size() > 0){
                        all_neighbors_are_dead_ends = false;
                    }
                }
                for (auto& neighbor : neighbors){
                    //choose the path with highest coverage
                    if (coverage[neighbor.first] < max_coverage/5.0
                            && coverage[neighbor.first] < coverage[contig_name]/5.0
                            && length_of_contigs[neighbor.first] < 2*length_of_contigs[contig_name]){
                        links_to_delete.insert({neighbor, {contig_name,end}});
                        links_to_delete.insert({{contig_name,end}, neighbor});
                    }
                    //if the neighbor is a dead end that looks like an error tip (short and less covered than both the best
                    //neighbor and the contig), cut it. A long or well-covered dead end is more likely a true contig whose
                    //other side is disconnected (e.g. a coverage gap at this k)
                    bool looks_like_error_tip = length_of_contigs[neighbor.first] < 3*k
                            && coverage[neighbor.first] < 0.5*std::min(max_coverage, coverage[contig_name]);
                    if (!all_neighbors_are_dead_ends && looks_like_error_tip &&
                            (links[neighbor.first][0].size() == 0 || links[neighbor.first][1].size() == 0)){
                        links_to_delete.insert({neighbor, {contig_name,end}});
                        links_to_delete.insert({{contig_name,end}, neighbor});
                    }
                }
            }
        }
        
    }

    //now write the gfa file without the links to delete
    ifstream input(gfa_in);
    string line;
    ofstream out(gfa_out);
    while (std::getline(input, line))
    {
        if (line[0] == 'L'){
            string name1;
            string name2;
            string orientation1;
            string orientation2;
            string dont_care;
            std::stringstream ss(line);
            int i = 0;
            ss >> dont_care >> name1 >> orientation1 >> name2 >> orientation2;

            if (links_to_delete.find({{graph.find(name1), (orientation1 == "+" ? 1 : 0)}, {graph.find(name2), (orientation2 == "+" ? 0 : 1)}}) == links_to_delete.end()){
                out << line << "\n";
            }
        }
        else{
            out << line << "\n";
        }
    }
}

/**
 * @brief simple function to trim tips and isolated nodes with a coverage below min_coverage and a length below min_length
 * 
 * @param gfa_in 
 * @param min_coverage 
 * @param min_length 
 * @param gfa_out 
 */
void trim_tips_isolated_contigs_and_bubbles(std::string gfa_in, int min_coverage, int min_length, std::string gfa_out, bool single_genome, bool hard_contiguity){
    LinkedGraph graph = load_linked_graph(gfa_in, false);
    check_depths(graph);
    check_links(graph);
    const vector<std::array<vector<pair<int,char>>, 2>>& links = graph.links;
    const vector<float>& coverage = graph.depth;
    const vector<int>& length_of_contigs = graph.length;

    //now trim the tips, isolated contigs and bubbles with a coverage below min_coverage and a length below min_length (bubbles need to have coverage 1)
    vector<bool> contigs_to_remove (graph.size(), false);
    //links of contigs that look valid but probably reduce contiguity, to detach from the graph: ((contig, end), (neighbor, end of neighbor))
    std::set<pair<pair<int,char>, pair<int,char>>> links_to_detach;
    for (int contig_name = 0 ; contig_name < graph.size() ; contig_name++){
        //remove tip or isolated contig
        if (links[contig_name][0].size() == 0 || links[contig_name][1].size() == 0){
            if (coverage[contig_name] < min_coverage && (length_of_contigs[contig_name] < min_length || coverage[contig_name]==1)){
                contigs_to_remove[contig_name] = true;
            }
            else { //the contig is valid but detaching it may improve the contiguity
                for (const auto& l : links[contig_name][0]) {
                    if (coverage[l.first] > 2 * coverage[contig_name]) {
                        links_to_detach.insert({{contig_name, 0}, l});
                    }
                }
                for (const auto& l : links[contig_name][1]) {
                    if (coverage[l.first] > 2 * coverage[contig_name]) {
                        links_to_detach.insert({{contig_name, 1}, l});
                    }
                }
            }
        }
        //remove bubble
        if (links[contig_name][0].size() == 1 && links[contig_name][1].size() == 1) {
            int neighbor_left = links[contig_name][0][0].first;
            char end_of_neighbor_left = links[contig_name][0][0].second;
            int neighbor_right = links[contig_name][1][0].first;
            char end_of_neighbor_right = links[contig_name][1][0].second;


            if (coverage[contig_name] == 1 || hard_contiguity) {
                int other_neighbor_of_contig_left = -1;
                if (end_of_neighbor_left == 0 && links[neighbor_left][0].size() == 2) {
                    for (const auto& l : links[neighbor_left][0]) {
                        if (l.first != contig_name) {
                            other_neighbor_of_contig_left = l.first;
                        }
                    }
                } else if (end_of_neighbor_left == 1 && links[neighbor_left][1].size() == 2) {
                    for (const auto& l : links[neighbor_left][1]) {
                        if (l.first != contig_name) {
                            other_neighbor_of_contig_left = l.first;
                        }
                    }
                }

                int other_neighbor_of_contig_right = -1;
                if (end_of_neighbor_right == 0 && links[neighbor_right][0].size() == 2) {
                    for (const auto& l : links[neighbor_right][0]) {
                        if (l.first != contig_name) {
                            other_neighbor_of_contig_right = l.first;
                        }
                    }
                } else if (end_of_neighbor_right == 1 && links[neighbor_right][1].size() == 2) {
                    for (const auto& l : links[neighbor_right][1]) {
                        if (l.first != contig_name) {
                            other_neighbor_of_contig_right = l.first;
                        }
                    }
                }

                if (other_neighbor_of_contig_left == other_neighbor_of_contig_right && other_neighbor_of_contig_left != -1) {
                    // This is a real bubble with two paths between same pair of nodes
                    if (coverage[contig_name] < coverage[other_neighbor_of_contig_left] && coverage[contig_name] < min_coverage) {
                        contigs_to_remove[contig_name] = true;
                    }
                    else if (hard_contiguity && coverage[contig_name]*length_of_contigs[contig_name] < coverage[other_neighbor_of_contig_left]*length_of_contigs[other_neighbor_of_contig_left]){
                        links_to_detach.insert({{contig_name, 0}, links[contig_name][0][0]});
                        links_to_detach.insert({{contig_name, 1}, links[contig_name][1][0]});
                    }
                }
            }
        }
        //if contig attached among contigs of higher coverage, detach
        if (links[contig_name][0].size() > 0 && links[contig_name][1].size() > 0) {
            // Check if all neighbors have significantly higher coverage
            bool all_neighbors_higher_coverage = true;
            for (const auto& l : links[contig_name][0]) {
                if (coverage[l.first] < 2 * coverage[contig_name]) {
                    all_neighbors_higher_coverage = false;
                    break;
                }
            }
            if (all_neighbors_higher_coverage) {
                for (const auto& l : links[contig_name][1]) {
                    if (coverage[l.first] < 2 * coverage[contig_name]) {
                        all_neighbors_higher_coverage = false;
                        break;
                    }
                }
            }

            // If all neighbors have higher coverage, check if they have alternative neighbors better covered
            if (all_neighbors_higher_coverage) {
                bool all_neighbors_have_alternatives = true;
                
                // Check left neighbors
                for (const auto& l : links[contig_name][0]) {
                    int neighbor = l.first;
                    char neighbor_end = l.second;
                    auto& neighbor_links = (neighbor_end == 0) ? links[neighbor][0] : links[neighbor][1];
                    
                    bool has_better_alternative = false;
                    for (const auto& alt_link : neighbor_links) {
                        if (alt_link.first != contig_name && coverage[alt_link.first] >= 2*coverage[contig_name]) {
                            has_better_alternative = true;
                            break;
                        }
                    }
                    if (!has_better_alternative) {
                        all_neighbors_have_alternatives = false;
                        break;
                    }
                }
                
                // Check right neighbors
                if (all_neighbors_have_alternatives) {
                    for (const auto& l : links[contig_name][1]) {
                        int neighbor = l.first;
                        char neighbor_end = l.second;
                        auto& neighbor_links = (neighbor_end == 0) ? links[neighbor][0] : links[neighbor][1];
                        
                        bool has_better_alternative = false;
                        for (const auto& alt_link : neighbor_links) {
                            if (alt_link.first != contig_name && coverage[alt_link.first] >= 2*coverage[contig_name]) {
                                has_better_alternative = true;
                                break;
                            }
                        }
                        if (!has_better_alternative) {
                            all_neighbors_have_alternatives = false;
                            break;
                        }
                    }
                }
                
                if (all_neighbors_have_alternatives) {
                    // Detach all links of the contig
                    for (const auto& l : links[contig_name][0]) {
                        links_to_detach.insert({{contig_name, 0}, l});
                    }
                    for (const auto& l : links[contig_name][1]) {
                        links_to_detach.insert({{contig_name, 1}, l});
                    }
                }
            }
        }
    }

    //now write the gfa file without the contigs to remove
    ifstream input(gfa_in);
    string line;
    ofstream out(gfa_out);
    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            if (!contigs_to_remove[graph.find(name)]){
                out << line << "\n";
            }
        }
        else if (line[0] == 'L'){
            string name1;
            string name2;
            string orientation1;
            string orientation2;
            string dont_care;
            std::stringstream ss(line);
            int i = 0;
            ss >> dont_care >> name1 >> orientation1 >> name2 >> orientation2;

            pair<int,char> end1 = {graph.find(name1), (orientation1 == "+" ? 1 : 0)};
            pair<int,char> end2 = {graph.find(name2), (orientation2 == "+" ? 0 : 1)};
            if (!contigs_to_remove[graph.find(name1)]
                && !contigs_to_remove[graph.find(name2)]
                && links_to_detach.find({end1, end2}) == links_to_detach.end() 
                && links_to_detach.find({end2, end1}) == links_to_detach.end()){
                out << line << "\n";
            }
        }
        else{
            out << line << "\n";
        }
    }

}

void load_GFA(string gfa_file, vector<Segment> &segments, unordered_map<string, int> &segment_IDs, bool load_in_RAM){
    //load the segments from the GFA file
    
    //in a first pass index all the segments by their name
    ifstream gfa(gfa_file);
    string line;
    long int next_pos_in_file = 0;
    while (getline(gfa, line)){
        long int pos_in_file = next_pos_in_file;
        next_pos_in_file += line.size() + 1;
        if (line[0] == 'S'){
            stringstream ss(line);
            string nothing, name, seq;
            ss >> nothing >> name >> seq;

            double coverage = 0;
            //try to find a DP: tag
            string tag;
            while (ss >> tag){
                if (tag.substr(0, 3) == "DP:" || tag.substr(0, 3) == "dp:" || tag.substr(0, 3) == "km:"){
                    coverage = std::stof(tag.substr(5, tag.size() - 5));
                }
            }

            if (load_in_RAM == false){
                Segment s(name, segments.size(), vector<pair<vector<pair<int,int>>, vector<string>>>(2), pos_in_file, seq.size(), coverage);
                segment_IDs[name] = s.ID;
                segments.push_back(std::move(s));
            }
            else{
                Segment s(name, segments.size(), vector<pair<vector<pair<int,int>>, vector<string>>>(2), pos_in_file, seq, seq.size(), coverage);
                segment_IDs[name] = s.ID;
                segments.push_back(std::move(s));
            }
        }
    }
    gfa.close();

    //in a second pass, load the links
    gfa.open(gfa_file);
    while (getline(gfa, line)){
        if (line[0] == 'L'){
            stringstream ss(line);
            string nothing, name1, name2;
            string orientation1, orientation2;
            string cigar;

            ss >> nothing >> name1 >> orientation1 >> name2 >> orientation2 >> cigar;

            int end1 = 1;
            int end2 = 0;

            if (orientation1 == "-"){
                end1 = 0;
            }
            if (orientation2 == "-"){
                end2 = 1;
            }

            auto it1 = segment_IDs.find(name1);
            auto it2 = segment_IDs.find(name2);
            if (it1 == segment_IDs.end() || it2 == segment_IDs.end()){
                cerr << "WARNING: link between unknown segments ignored in " << gfa_file << ": " << line << "\n";
                continue;
            }
            int ID1 = it1->second;
            int ID2 = it2->second;

            //check that the link did not already exist
            bool already_exists = false;
            for (const pair<int,int>& link : segments[ID1].links[end1].first){
                if (link.first == ID2 && link.second == end2){
                    already_exists = true;
                }
            }
            if (!already_exists){

                segments[ID1].links[end1].first.push_back({ID2, end2});
                segments[ID1].links[end1].second.push_back(cigar);

                segments[ID2].links[end2].first.push_back({ID1, end1});
                segments[ID2].links[end2].second.push_back(cigar);
            }
        }
    }
    gfa.close();
}

/**
 * @brief Merge all old_segments into new_segments
 * 
 * @param old_segments 
 * @param new_segments 
 * @param original_gfa_file original gfa file to retrieve the sequences
 * @param rename rename the contigs in short names or keep the old names with underscores in between
 * @return * void 
 */
void merge_adjacent_contigs(vector<Segment> &old_segments, vector<Segment> &new_segments, string original_gfa_file, bool rename, int num_threads){

    //old IDs of segments that have already been looked at and merged (don't want to merge them twice). Atomic because it is read outside of the critical section
    vector<std::atomic<bool>> already_looked_at_segments (old_segments.size());
    for (auto& a : already_looked_at_segments){
        a = false;
    }
    unordered_map<pair<int,int>,pair<int,int>, PairHash> old_ID_to_new_ID; //associates (old_id, old end) with (new_id, new_end)
    int number_of_merged_contigs = 0;
    set<pair<pair<pair<int,int>, pair<int,int>>,string>> links_to_add; //list of links to add, all in old IDs and old ends
    omp_lock_t lock_new_segment; //locks the creating of new segments, including the additions to links_to_add
    omp_init_lock(&lock_new_segment);

    double total_time_pre = 0;
    double total_time_prepare = 0;
    double total_time_merge = 0;

    vector<vector<string>> all_names (num_threads);
    vector<vector<string>> all_seqs (num_threads);
    vector<vector<double>> all_coverages (num_threads);
    vector<vector<int>> all_lengths(num_threads);
    vector<vector<int>> all_IDs(num_threads);
    for (int t = 0 ; t < num_threads ; t++){
        all_names[t].reserve(100);
        all_seqs[t].reserve(100);
        all_coverages[t].reserve(100);
        all_lengths[t].reserve(100);
        all_IDs[t].reserve(100);
    }

    #pragma omp parallel for num_threads(num_threads)
    for (int seg_idx = 0 ; seg_idx < old_segments.size() ; seg_idx++){

        int thread_num = omp_get_thread_num();

        auto time_start = std::chrono::high_resolution_clock::now();
        Segment& old_seg = old_segments[seg_idx];

        if (already_looked_at_segments[old_seg.ID]){
            continue;
        }
        //check if it has either at least two neighbors left or that its neighbor left has at least two neighbors right
        bool dead_end_left = false;
        if (old_seg.links[0].first.size() != 1 || old_segments[old_seg.links[0].first[0].first].links[old_seg.links[0].first[0].second].first.size() != 1 || old_segments[old_seg.links[0].first[0].first].ID == old_seg.ID){
            dead_end_left = true;
        }

        bool dead_end_right = false;
        if (old_seg.links[1].first.size() != 1 || old_segments[old_seg.links[1].first[0].first].links[old_seg.links[1].first[0].second].first.size() != 1 || old_segments[old_seg.links[1].first[0].first].ID == old_seg.ID){
            dead_end_right = true;
        }

        if (!dead_end_left && !dead_end_right){ //means this contig is in the middle of a long haploid contig, no need to merge
            continue;
        }

        auto time_before_prepare = std::chrono::high_resolution_clock::now();

        int other_end_of_merged_contig_ID = old_seg.ID;
        int other_end_of_merged_contig_end = 0;

        //prepare the merge of the contig (don't merge yet to avoid conflicts with other threads)
        all_IDs[thread_num] = {old_seg.ID};
        all_names[thread_num].clear();
        all_seqs[thread_num].clear();
        all_coverages[thread_num].clear();
        all_lengths[thread_num].clear();
        if (dead_end_left && !dead_end_right){
            //let's see how far we can go right
            all_names[thread_num] = {old_seg.name};
            all_seqs[thread_num] = {old_seg.get_seq(original_gfa_file)};
            all_coverages[thread_num] = {old_seg.get_coverage()};
            all_lengths[thread_num] = {old_seg.get_length()};
            int current_ID = old_seg.ID;
            int current_end = 1;

            while (old_segments[current_ID].links[current_end].first.size() == 1 && old_segments[old_segments[current_ID].links[current_end].first[0].first].links[old_segments[current_ID].links[current_end].first[0].second].first.size() == 1){
                string cigar = old_segments[current_ID].links[current_end].second[0];
                int tmp_current_end = 1-old_segments[current_ID].links[current_end].first[0].second;
                current_ID = old_segments[current_ID].links[current_end].first[0].first;
                current_end = tmp_current_end;
                all_names[thread_num].push_back(old_segments[current_ID].name);
                string seq = old_segments[current_ID].get_seq(original_gfa_file);

                //now reverse complement if current_end is 1
                if (current_end == 0){
                    seq = reverse_complement(seq);
                }
                //trim the sequence if there is a CIGAR
                size_t num_matches = min((size_t) overlap_length_of_CIGAR(cigar), seq.size());
                all_seqs[thread_num].push_back(seq.substr(num_matches));
                all_coverages[thread_num].push_back(old_segments[current_ID].get_coverage());
                all_lengths[thread_num].push_back(old_segments[current_ID].get_length());
                all_IDs[thread_num].push_back(current_ID);
            }

            other_end_of_merged_contig_ID = current_ID;
            other_end_of_merged_contig_end = current_end;
        }
        else if (!dead_end_left && dead_end_right){         
            //let's see how far we can go left
            all_names[thread_num] = {old_seg.name};
            string seq = old_seg.get_seq(original_gfa_file);
            all_seqs[thread_num] = {reverse_complement(seq)};
            all_coverages[thread_num] = {old_seg.get_coverage()};
            all_lengths[thread_num] = {old_seg.get_length()};
            int current_ID = old_seg.ID;
            int current_end = 0;

            // cout << "exploring all the contigs left" << endl;
            // cout << "first exploring the link between " << old_segments[current_ID].name << " and " << old_segments[old_segments[current_ID].links[current_end].first[0].first].name << endl;
            
            while (old_segments[current_ID].links[current_end].first.size() == 1 && old_segments[old_segments[current_ID].links[current_end].first[0].first].links[old_segments[current_ID].links[current_end].first[0].second].first.size() == 1){
                string cigar = old_segments[current_ID].links[current_end].second[0];
                int tmp_current_end = 1-old_segments[current_ID].links[current_end].first[0].second;
                current_ID = old_segments[current_ID].links[current_end].first[0].first;
                current_end = tmp_current_end;
                all_names[thread_num].push_back(old_segments[current_ID].name);
                string seq = old_segments[current_ID].get_seq(original_gfa_file);
                //now reverse complement if current_end is 0
                if (current_end == 0){
                    seq = reverse_complement(seq);
                }
                //trim the sequence if there is a CIGAR
                size_t num_matches = min((size_t) overlap_length_of_CIGAR(cigar), seq.size());
                all_seqs[thread_num].push_back(seq.substr(num_matches));
                all_coverages[thread_num].push_back(old_segments[current_ID].get_coverage());
                all_lengths[thread_num].push_back(old_segments[current_ID].get_length());
                all_IDs[thread_num].push_back(current_ID);
            }

            other_end_of_merged_contig_ID = current_ID;
            other_end_of_merged_contig_end = current_end;
        }
        
        auto time_after_prepare = std::chrono::high_resolution_clock::now();

        //check that we can proceed thread-safely
        bool thread_safe = true;
        #pragma omp critical
        {
            if (already_looked_at_segments[old_seg.ID] || already_looked_at_segments[other_end_of_merged_contig_ID]){
                thread_safe = false;
            }
            else{
                for (int ID : all_IDs[thread_num]){
                    already_looked_at_segments[ID] = true;
                }
            }
        }

        //actually merge the contigs
        if (thread_safe){    
            if (dead_end_left && dead_end_right){

                omp_set_lock(&lock_new_segment);
                string name = old_seg.name;
                if (rename){
                    name = std::to_string(number_of_merged_contigs);
                    number_of_merged_contigs++;
                }

                new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), old_seg.get_length(), old_seg.get_coverage()));
                old_ID_to_new_ID[{old_seg.ID, 0}] = {new_segments.size() - 1, 0};
                old_ID_to_new_ID[{old_seg.ID, 1}] = {new_segments.size() - 1, 1};
                new_segments[new_segments.size()-1].seq = old_seg.get_seq(original_gfa_file);
            
                //add the links
                int idx_link = 0;
                for (pair<int,int> link : old_seg.links[0].first){
                    links_to_add.insert({{{old_seg.ID, 0}, link}, old_seg.links[0].second[idx_link]});
                    idx_link++;
                }
                idx_link = 0;
                for (pair<int,int> link : old_seg.links[1].first){
                    links_to_add.insert({{{old_seg.ID, 1}, link}, old_seg.links[1].second[idx_link]});
                    idx_link++;
                }
                omp_unset_lock(&lock_new_segment);
                
            }
            else if (dead_end_left && !dead_end_right){
                //create the new contig
                string new_name = "";
                for (const string& name : all_names[thread_num]){
                    new_name += name + "_";
                }
                new_name = new_name.substr(0, new_name.size()-1);
                string new_seq = "";
                for (const string& seq : all_seqs[thread_num]){
                    new_seq += seq;
                }
                double new_coverage = 0;
                int new_length = 0;
                int idx = 0;
                for (double coverage : all_coverages[thread_num]){
                    new_coverage += all_coverages[thread_num][idx]*all_lengths[thread_num][idx];
                    new_length += all_lengths[thread_num][idx];
                    idx++;
                }
                new_coverage = new_coverage/new_length;

                omp_set_lock(&lock_new_segment);
                string name = new_name;
                if (rename){
                    name = std::to_string(number_of_merged_contigs);
                    number_of_merged_contigs++;
                }

                new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), (int) new_seq.size(), new_coverage)); //coverage is weighted by the lengths of the original segments
                new_segments[new_segments.size()-1].seq = new_seq;
                old_ID_to_new_ID[{old_seg.ID, 0}] = {new_segments.size() - 1, 0};
                old_ID_to_new_ID[{other_end_of_merged_contig_ID, other_end_of_merged_contig_end}] = {new_segments.size() - 1, 1};

                //add the links
                int idx_link = 0;
                for (pair<int,int> link : old_seg.links[0].first){
                    links_to_add.insert({{{old_seg.ID, 0}, link}, old_seg.links[0].second[idx_link]});
                    idx_link++;
                }
                idx_link = 0;
                for (pair<int,int> link : old_segments[other_end_of_merged_contig_ID].links[other_end_of_merged_contig_end].first){
                    links_to_add.insert({{{other_end_of_merged_contig_ID, other_end_of_merged_contig_end}, link}, old_segments[other_end_of_merged_contig_ID].links[other_end_of_merged_contig_end].second[idx_link]});
                    idx_link++;
                }
                omp_unset_lock(&lock_new_segment);
            }
            else if (!dead_end_left && dead_end_right){
            
                //create the new contig
                string new_name = "r";
                for (const string& name : all_names[thread_num]){
                    new_name += name + "_";
                }
                new_name = new_name.substr(0, new_name.size()-1);
                string new_seq = "";
                for (const string& seq : all_seqs[thread_num]){
                    new_seq += seq;
                }
                double new_coverage = 0;
                int new_length = 0;
                int idx = 0;
                for (double coverage : all_coverages[thread_num]){
                    new_coverage += all_coverages[thread_num][idx]*all_lengths[thread_num][idx];
                    new_length += all_lengths[thread_num][idx];
                    idx++;
                }
                new_coverage = new_coverage/new_length;

                omp_set_lock(&lock_new_segment);  
                string name = new_name;
                if (rename){
                    name = std::to_string(number_of_merged_contigs);
                    number_of_merged_contigs++;
                }

                new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), (int) new_seq.size(), new_coverage)); //coverage is weighted by the lengths of the original segments
                new_segments[new_segments.size()-1].seq = new_seq;
                old_ID_to_new_ID[{old_seg.ID, 1}] = {new_segments.size() - 1, 0};
                old_ID_to_new_ID[{other_end_of_merged_contig_ID, other_end_of_merged_contig_end}] = {new_segments.size() - 1, 1};

                // add the links
                int idx_link = 0;
                for (pair<int,int> link : old_seg.links[1].first){
                    links_to_add.insert({{{old_seg.ID, 1}, link}, old_seg.links[1].second[idx_link]});
                    idx_link++;
                }
                idx_link = 0;
                for (pair<int,int> link : old_segments[other_end_of_merged_contig_ID].links[other_end_of_merged_contig_end].first){
                    links_to_add.insert({{{other_end_of_merged_contig_ID, other_end_of_merged_contig_end}, link}, old_segments[other_end_of_merged_contig_ID].links[other_end_of_merged_contig_end].second[idx_link]});
                    idx_link++;
                }
                omp_unset_lock(&lock_new_segment);
            }
        }

        auto time_after_merge = std::chrono::high_resolution_clock::now();

        total_time_pre = std::chrono::duration_cast<std::chrono::milliseconds>(time_before_prepare-time_start).count();
        total_time_prepare = std::chrono::duration_cast<std::chrono::milliseconds>(time_after_prepare - time_before_prepare).count();
        total_time_merge = std::chrono::duration_cast<std::chrono::milliseconds>(time_after_merge - time_after_prepare).count();
    }
    omp_destroy_lock(&lock_new_segment);

    //some contigs are left: the ones that were in circular rings... go through them and add them
    for (Segment& old_seg : old_segments){
        if (!already_looked_at_segments[old_seg.ID]){
            int current_ID = old_seg.ID;
            int current_end = 1;
            vector<string> all_names = {old_seg.name};
            vector<string> all_seqs = {old_seg.get_seq(original_gfa_file)};
            vector<double> all_coverages = {old_seg.get_coverage()};
            vector<int> all_lengths = {old_seg.get_length()};
            bool circular_as_expected = true;
            while (old_segments[current_ID].links[current_end].first.size() == 1 && old_segments[old_segments[current_ID].links[current_end].first[0].first].links[old_segments[current_ID].links[current_end].first[0].second].first.size() == 1){
                if (old_segments[current_ID].links[current_end].first[0].first == old_seg.ID){
                    break;
                }
                already_looked_at_segments[current_ID] = true;
                string cigar = old_segments[current_ID].links[current_end].second[0];
                int tmp_current_end = 1-old_segments[current_ID].links[current_end].first[0].second;
                current_ID = old_segments[current_ID].links[current_end].first[0].first;
                current_end = tmp_current_end;
                all_names.push_back(old_segments[current_ID].name);
                string seq = old_segments[current_ID].get_seq(original_gfa_file);
                //now reverse complement if current_end is 0
                if (current_end == 0){
                    seq = reverse_complement(seq);
                }
                //trim the sequence if there is a CIGAR
                size_t num_matches = min((size_t) overlap_length_of_CIGAR(cigar), seq.size());
                all_seqs.push_back(seq.substr(num_matches));
                all_coverages.push_back(old_segments[current_ID].get_coverage());
                all_lengths.push_back(old_segments[current_ID].get_length());
            }
            if (old_segments[current_ID].links[current_end].first.size() != 1 || old_segments[current_ID].links[current_end].first[0].first != old_seg.ID){
                circular_as_expected = false;
            }
            if (circular_as_expected){
                already_looked_at_segments[current_ID] = true;
                string new_name = "";
                for (const string& name : all_names){
                    new_name += name + "_";
                }
                new_name = new_name.substr(0, new_name.size()-1);
                string new_seq = "";
                for (const string& seq : all_seqs){
                    new_seq += seq;
                }
                double new_coverage = 0;
                int new_length = 0;
                int idx = 0;
                for (double coverage : all_coverages){
                    new_coverage += all_coverages[idx]*all_lengths[idx];
                    new_length += all_lengths[idx];
                    idx++;
                }
                new_coverage = new_coverage/new_length;
                string name = new_name;
                if (rename){

                    name = std::to_string(number_of_merged_contigs);
                    number_of_merged_contigs++;
                }
                new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), (int) new_seq.size(), new_coverage)); //coverage is weighted by the lengths of the original segments
                new_segments[new_segments.size()-1].seq = new_seq;

                old_ID_to_new_ID[{old_seg.ID, 0}] = {new_segments.size() - 1, 0};
                old_ID_to_new_ID[{current_ID, current_end}] = {new_segments.size() - 1, 1};

                //add the link to circularize
                int idx_link = 0;
                for (pair<int,int> link : old_seg.links[0].first){
                    links_to_add.insert({{{old_seg.ID, 0}, link}, old_seg.links[0].second[idx_link]});
                    idx_link++;
                }
            }
            else{
                cout << "ERROR: contig " << old_seg.name << " was discarded while merging the reads in graphunzip.cpp" << endl;
            }
        }
    }

    //now add the links in the new segments (skipping the links to segments that were discarded)
    for (const pair<pair<pair<int,int>, pair<int,int>>, string>& link : links_to_add){
        auto from = old_ID_to_new_ID.find(link.first.first);
        auto to = old_ID_to_new_ID.find(link.first.second);
        if (from == old_ID_to_new_ID.end() || to == old_ID_to_new_ID.end()){
            continue;
        }
        new_segments[from->second.first].links[from->second.second].first.push_back(to->second);
        new_segments[from->second.first].links[from->second.second].second.push_back(link.second);
    }
}

void output_graph(string gfa_output, string gfa_input, vector<Segment> &segments){
    ofstream gfa(gfa_output);
    for (Segment& s : segments){
        if (s.name != "delete_me"){
            gfa << "S\t" << s.name << "\t" << s.get_seq(gfa_input) << "\tDP:f:" << s.get_coverage() <<  "\n";
        }
    }
    for (Segment& s : segments){
        for (int end = 0 ; end < 2 ; end++){
            for (int neigh = 0 ; neigh < s.links[end].first.size() ; neigh++){

                //to make sure the link is not outputted twice
                if (s.ID > s.links[end].first[neigh].first || (s.ID == s.links[end].first[neigh].first && end > s.links[end].first[neigh].second) ){
                    continue;
                }
                if (s.name == "delete_me" || segments[s.links[end].first[neigh].first].name == "delete_me"){
                    continue;
                }

                string orientation = "+";
                if (end == 0){
                    orientation = "-";
                }
                gfa << "L\t" << s.name << "\t" << orientation << "\t" << segments[s.links[end].first[neigh].first].name << "\t";
                if (s.links[end].first[neigh].second == 0){
                    gfa << "+\t";
                }
                else{
                    gfa << "-\t";
                }
                gfa << s.links[end].second[neigh] << "\n";
            }
        }
    }
    gfa.close();
}

