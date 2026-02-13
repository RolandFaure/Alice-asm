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

void shave(std::string input_file, std::string output_file, int max_length){
    std::ifstream input(input_file);
    if (!input.is_open())
    {
        std::cout << "Could not open file iicy " << input_file << std::endl;
        exit(1);
    }

    std::string line;
    set<string> good_contigs;
    unordered_map<string, pair<bool,bool>> linked;

    std::ofstream output(output_file);

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

            if (linked.find(name1) == linked.end()){
                linked[name1] = {false, false};
            }
            if (linked.find(name2) == linked.end()){
                linked[name2] = {false, false};
            }
            if (orientation1 == "+"){
                linked[name1].second = true;
            }
            else{
                linked[name1].first = true;
            }

            if (orientation2 == "+"){
                linked[name2].first = true;
            }
            else{
                linked[name2].second = true;
            }
        }
    }

    for (auto i: linked){
        if (i.second.first && i.second.second){
            good_contigs.insert(i.first);
        }
    }

    input.close();
    input.open(input_file);

    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            if (sequence.size() > max_length || good_contigs.find(name) != good_contigs.end()){
                output << line << "\n";
                good_contigs.insert(name);
            }
        }
        else if (line[0] == 'L'){
            string name1;
            string name2;
            string dont_care;
            std::stringstream ss(line);
            int i = 0;
            ss >> dont_care >> name1 >> dont_care >> name2;
            // std::cerr << name1 << " " << name2 << "\n";
            // for (auto i: good_contigs){
            //     std::cerr << i << ",";
            // }
            // std::cerr << "\n";
            if (good_contigs.find(name1) != good_contigs.end() && good_contigs.find(name2) != good_contigs.end()){
                output << line << "\n";
            }
        }
        else{
            output << line << "\n";
        }
    }
}

void compute_exact_CIGARs(std::string gfa_in, std::string gfa_out, int max_overlap, int default_overlap){

    //go through the graph and for all links, compute the exact CIGAR (that will be only M)
    ifstream input(gfa_in);
    ofstream out(gfa_out);

    string line;
    //first index the position of every contig in the file
    unordered_map<string, long int> pos_of_contig_seq_in_file;
    long int pos = 0;
    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            pos_of_contig_seq_in_file[name] = pos;
        }
        pos += line.size() + 1;
    }

    input.close();
    input.open(gfa_in);
    ifstream input2(gfa_in);

    //now, for each L line, compute the exact CIGAR
    while (std::getline(input, line))
    {
        if (line[0] == 'L')
        {
            string name1;
            string name2;
            string orientation1;
            string orientation2;
            string dont_care;
            std::stringstream ss(line);
            int i = 0;
            ss >> dont_care >> name1 >> orientation1 >> name2 >> orientation2;

            //get the sequences of the two contigs
            input2.seekg(pos_of_contig_seq_in_file[name1]);
            std::getline(input2, line);
            string seq1, seq2;
            std::stringstream ss2(line);
            ss2 >> dont_care >> dont_care >> seq1;
            input2.seekg(pos_of_contig_seq_in_file[name2]);
            std::getline(input2, line);

            ss2 = std::stringstream(line);
            ss2 >> dont_care >> dont_care >> seq2;

            if (orientation1 == "-"){
                seq1 = reverse_complement(seq1);
            }
            if (orientation2 == "-"){
                seq2 = reverse_complement(seq2);
            }

            int overlap = 30;
            if (seq1.size() < 30 || seq2.size() < 30){ //can happen if they could not be reconstructed
                overlap = 0;
            }
            else{
                string end1 = seq1.substr(seq1.size()-30, 30); 
                while (overlap < seq1.size() && overlap < seq2.size() && overlap < max_overlap && end1 != seq2.substr(overlap-30, 30) ){
                    overlap++;
                    // cout << "overlapping " << end1 << " and " << seq2.substr(overlap-30, 30) << "\n";
                }
                if (overlap == seq1.size() || overlap == seq2.size() || overlap == max_overlap){
                    // cerr << "ERROR: no overlap found between " << name1 << " and " << name2 << "\n";
                    // exit(1);
                    overlap = default_overlap;
                }
            }
            // cout << "overlap between " << name1 << " and " << name2 << " is " << overlap << "\n";
            out << "L\t" << name1 << "\t" << orientation1 << "\t" << name2 << "\t" << orientation2 << "\t" << overlap << "M\n";

        }
        else{
            out << line << "\n";
        }
    }
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

    for (auto c: contigs){
        output << c << "\n";
    }
    for (auto l: links){
        output << l << "\n";
    }
    output.close();

    //move the sorted file to the original file
    std::string command = "mv " + gfa + ".sorted " + gfa;
    auto res = system(command.c_str());
}


struct Path{
    vector<string> contigs;
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
    
    unordered_flat_map<uint64_t, pair<string,int>> kmers_to_contigs; //in what contig is the kmer and at what position (only unique kmer ofc, meant to work with unitig graph)
    unordered_flat_map<string, int> length_of_contigs;
    unordered_flat_map<string, string> contig_sequences; // Store contig sequences in memory
    unordered_map<string, pair<vector<pair<string, char>>, vector<pair<string,char>>>> linked;

    ifstream input(unitig_graph);
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
            length_of_contigs[name] = sequence.size();
            contig_sequences[name] = sequence; // Store sequence in memory

            if (linked.find(name) == linked.end()){
                linked[name] = {vector<pair<string,char>>(0), vector<pair<string,char>>(0)};
            }

            uint64_t hash_foward = -1;
            size_t pos_end = 0;
            long pos_begin = -km;
            while (roll_f(hash_foward, km, sequence, pos_end, pos_begin, false)){
                
                if (pos_begin>=0  && pos_begin % 5 == 0){
                    kmers_to_contigs[hash_foward] = make_pair(name, pos_end-km);
                }
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

            char or1 = (orientation1 == "+" ? 1 : 0);
            char or2 = (orientation2 == "+" ? 0 : 1);
            auto neighbor = make_pair(name2, or2);
            if (orientation1 == "+"){
                if (std::find(linked[name1].second.begin(), linked[name1].second.end(), neighbor) == linked[name1].second.end()){
                    linked[name1].second.push_back(neighbor);
                }
            }
            else{
                if (std::find(linked[name1].first.begin(), linked[name1].first.end(), neighbor) == linked[name1].first.end()){
                    linked[name1].first.push_back(neighbor);
                }
            }

            neighbor = make_pair(name1, or1);
            if (orientation2 == "+"){
                if (std::find(linked[name2].first.begin(), linked[name2].first.end(), neighbor) == linked[name2].first.end()){
                    linked[name2].first.push_back(neighbor);
                }
            }
            else{
                if (std::find(linked[name2].second.begin(), linked[name2].second.end(), neighbor) == linked[name2].second.end()){
                    linked[name2].second.push_back(neighbor);
                }
            }
        }
    }
    input.close();

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
    const int CHUNK_SIZE = 100;
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
            string previous_contig = "";
            bool previous_orientation = true;
            int previous_match_pos_on_contig = -1;
            bool read_is_single_path = true;
            
            while(roll(hash_foward, hash_reverse, km, line, pos_end, pos_begin, pos_middle, false)){
                if (pos_begin == pos_to_look_at){

                    unsigned long kmer = hash_foward; 
                    bool found_match = false;
                    bool forward_orientation = true;

                    if (kmers_to_contigs.find(kmer) != kmers_to_contigs.end()){
                        found_match = true;
                        forward_orientation = true;
                    }
                    else if (kmers_to_contigs.find(hash_reverse) != kmers_to_contigs.end()){
                        found_match = true;
                        forward_orientation = false;
                        kmer = hash_reverse;
                    }

                    if (found_match){
                        string contig = kmers_to_contigs[kmer].first;
                        int pos_in_contig = kmers_to_contigs[kmer].second;
                        
                        if (!forward_orientation){
                            pos_in_contig = length_of_contigs[contig] - pos_in_contig - km;
                        }

                        // Handle sequence before this match
                        if (previous_match_pos == -1){
                            // First match: add pos_begin bases of contig sequence
                            string contig_seq = contig_sequences[contig];
                            if (!forward_orientation) {
                                contig_seq = reverse_complement(contig_seq);
                            }
                            corrected_seq += contig_seq.substr(std::max((long int) 0, pos_in_contig - pos_begin), min(pos_begin, (long int) pos_in_contig));
                            if (pos_in_contig - pos_begin < 0 && !hard_correct){
                                corrected_seq = line.substr(0, pos_begin - pos_in_contig) + corrected_seq;
                            }
                            last_corrected_pos = pos_begin;
                        }
                        else if (contig == previous_contig && forward_orientation == previous_orientation){ //do a special case because that is very frequent
                            // Add sequence between previous match and current match on the same contig
                            string current_seq = contig_sequences[contig];
                            if (!forward_orientation) {
                                current_seq = reverse_complement(current_seq);
                            }
                            // Extract the portion between the two matches
                            corrected_seq += current_seq.substr(previous_match_pos_on_contig, pos_in_contig - previous_match_pos_on_contig);
                            last_corrected_pos = pos_begin;
                            // Update coverage for this portion
                            omp_set_lock(&coverage_lock);
                            coverages[contig] += (pos_in_contig - previous_match_pos_on_contig) / (double)length_of_contigs[contig];
                            omp_unset_lock(&coverage_lock);
                        }
                        else{
                            // Try to bridge from previous match to current match
                            int distance = max((long int) 1,pos_begin - previous_match_pos - pos_in_contig); // is distance is negative, still look the neighboring contig

                            vector<vector<pair<string, bool>>> all_paths;
                            if (distance > 0){
                                all_paths = list_all_paths_from_contig(linked, contig, !forward_orientation, distance, length_of_contigs, km);
                            }
                            else{
                                all_paths = {{{contig, !forward_orientation}}};
                            }
                            
                            int valid_path_count = 0;
                            int valid_path_index = -1;
                            for (int path_idx = 0; path_idx < all_paths.size(); path_idx++) {
                                auto node = all_paths[path_idx][all_paths[path_idx].size()-1];
                                if (node.first == previous_contig) {
                                    valid_path_count++;
                                    valid_path_index = path_idx;
                                }
                            }
                            
                            if (valid_path_count == 1) {
                                // Found unique path - use bridging contigs
                                const auto& bridging_path = all_paths[valid_path_index];
                                string correct_seq = "";

                                // First append to the read the end of the contig that was mathched last (previous_contig)
                                string prev_seq = contig_sequences[previous_contig]; 
                                if (!previous_orientation){
                                    prev_seq = reverse_complement(prev_seq);
                                }
                                correct_seq += prev_seq.substr(previous_match_pos_on_contig, prev_seq.size()-previous_match_pos_on_contig);
                                // Update coverage for the previous contig based on the portion used
                                omp_set_lock(&coverage_lock);
                                int length_used = min((long int) prev_seq.size()-previous_match_pos_on_contig, pos_begin-previous_match_pos+km-1);
                                coverages[previous_contig] += length_used / (double)length_of_contigs[previous_contig];
                                omp_unset_lock(&coverage_lock);

                                for (int i = bridging_path.size() - 2; i >= 1; i--) {
                                    string bridge_seq = contig_sequences[bridging_path[i].first];
                                    if (bridging_path[i].second) {
                                        bridge_seq = reverse_complement(bridge_seq);
                                    }
                                    // Skip overlap
                                    correct_seq += bridge_seq.substr(km-1);
                                    
                                    omp_set_lock(&coverage_lock);
                                    coverages[bridging_path[i].first] += 1;
                                    omp_unset_lock(&coverage_lock);
                                }

                                if (bridging_path.size() > 1){
                                    // Add the beginning of current contig to complete the bridging
                                    string current_seq = contig_sequences[contig];
                                    if (!forward_orientation) {
                                        current_seq = reverse_complement(current_seq);
                                    }
                                    // Add up to pos_in_contig bases from the current contig
                                    correct_seq += current_seq.substr(km-1, pos_in_contig);
                                    omp_set_lock(&coverage_lock);
                                    int length_used = pos_in_contig;
                                    coverages[contig] += length_used / (double)length_of_contigs[contig];
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

                        int length_of_contig_left = length_of_contigs[contig] - pos_in_contig - km;
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
                    string last_seq = contig_sequences[previous_contig];
                    int remaining_contig_length = length_of_contigs[previous_contig] - previous_match_pos_on_contig;
                    if (!previous_orientation){
                        last_seq = reverse_complement(last_seq);
                    }
                    int remaining_read_length = line.size() - last_corrected_pos;
                    corrected_seq += last_seq.substr(previous_match_pos_on_contig, min(remaining_contig_length, remaining_read_length));
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

    add_coverages_to_graph(unitig_graph, coverages);

    cout << "    -> Number of cleanly corrected/aligned reads: " << nb_reads_single_path << " out of " << nb_reads << endl;
}

void create_gaf_from_unitig_graph(std::string unitig_graph, int km, std::string reads_file, std::string output_file, robin_hood::unordered_flat_map<std::string, float>& coverages, int num_threads){
    
    unordered_flat_map<uint64_t, pair<string,int>> kmers_to_contigs;
    unordered_flat_map<string, int> length_of_contigs;
    unordered_map<string, pair<vector<pair<string, char>>, vector<pair<string,char>>>> linked;

    ifstream input(unitig_graph);
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
            length_of_contigs[name] = sequence.size();

            if (linked.find(name) == linked.end()){
                linked[name] = {vector<pair<string,char>>(0), vector<pair<string,char>>(0)};
            }

            uint64_t hash_foward = -1;
            size_t pos_end = 0;
            long pos_begin = -km;
            while (roll_f(hash_foward, km, sequence, pos_end, pos_begin, false)){
                
                if (pos_begin>=0  && pos_begin % 5 == 0){
                    kmers_to_contigs[hash_foward] = make_pair(name, pos_end-km);
                }
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

            char or1 = (orientation1 == "+" ? 1 : 0);
            char or2 = (orientation2 == "+" ? 0 : 1);
            auto neighbor = make_pair(name2, or2);
            if (orientation1 == "+"){
                if (std::find(linked[name1].second.begin(), linked[name1].second.end(), neighbor) == linked[name1].second.end()){
                    linked[name1].second.push_back(neighbor);
                }
            }
            else{
                if (std::find(linked[name1].first.begin(), linked[name1].first.end(), neighbor) == linked[name1].first.end()){
                    linked[name1].first.push_back(neighbor);
                }
            }

            neighbor = make_pair(name1, or1);
            if (orientation2 == "+"){
                if (std::find(linked[name2].first.begin(), linked[name2].first.end(), neighbor) == linked[name2].first.end()){
                    linked[name2].first.push_back(neighbor);
                }
            }
            else{
                if (std::find(linked[name2].second.begin(), linked[name2].second.end(), neighbor) == linked[name2].second.end()){
                    linked[name2].second.push_back(neighbor);
                }
            }
        }
    }
    input.close();

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
    const int CHUNK_SIZE = 100;
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
            string previous_contig = "";
            
            while(roll(hash_foward, hash_reverse, km, line, pos_end, pos_begin, pos_middle, false)){
                if (pos_begin == pos_to_look_at){

                    unsigned long kmer = hash_foward; 

                    if (kmers_to_contigs.find(kmer) != kmers_to_contigs.end()){
                        string contig = kmers_to_contigs[kmer].first;
                        int pos_in_contig = kmers_to_contigs[kmer].second;

                        if (current_path.contigs.size() == 0){
                            current_path.start_position_on_contig = pos_in_contig;
                        }
                        
                        if (current_path.contigs.size() > 0 && current_path.contigs[current_path.contigs.size()-1] == contig && current_path.orientations[current_path.orientations.size()-1] == true){
                            current_path.end_position_on_contig = pos_in_contig + km;
                        }
                        else{
                            if (previous_contig != ""){
                                int distance = pos_to_look_at - previous_match - pos_in_contig;
                                vector<vector<pair<string, bool>>> all_paths = list_all_paths_from_contig(linked, contig, false, distance, length_of_contigs, km);
                                
                                int valid_path_count = 0;
                                int valid_path_index = -1;
                                for (int path_idx = 0; path_idx < all_paths.size(); path_idx++) {
                                    auto node = all_paths[path_idx][all_paths[path_idx].size()-1];
                                    if (node.first == previous_contig) {
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
                                            coverages[path[i].first] += min(1.0, (line.size() - pos_begin) / (double)length_of_contigs[path[i].first]);
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
                            coverages[contig]+= std::min(1.0, (line.size()- pos_begin) / (double)length_of_contigs[contig]);
                            omp_unset_lock(&coverage_lock);
                        }
                        
                        int length_of_contig_left = length_of_contigs[contig] - pos_in_contig - km;
                        if (length_of_contig_left > 10){
                            pos_to_look_at += min((int) (length_of_contig_left*0.8) , (int)(line.size() - pos_begin - km - 5));
                        }
                        previous_match = pos_begin;
                        previous_contig = contig;
                    }
                    else if (kmers_to_contigs.find(hash_reverse) != kmers_to_contigs.end()){
                        string contig = kmers_to_contigs[hash_reverse].first;
                        int pos_in_contig = kmers_to_contigs[hash_reverse].second;

                        if (current_path.contigs.size() == 0){
                            current_path.start_position_on_contig = length_of_contigs[contig] - pos_in_contig - km;
                        }
                        
                        if (current_path.contigs.size() > 0 && current_path.contigs[current_path.contigs.size()-1] == contig && current_path.orientations[current_path.orientations.size()-1] == false){
                            current_path.end_position_on_contig = length_of_contigs[contig] - pos_in_contig;
                        }
                        else{
                            if (previous_contig != ""){
                                int distance = pos_to_look_at - previous_match - (length_of_contigs[contig] - pos_in_contig - km);
                                vector<vector<pair<string, bool>>> all_paths = list_all_paths_from_contig(linked, contig, true, distance, length_of_contigs, km);

                                int valid_path_count = 0;
                                int valid_path_index = -1;
                                for (int path_idx = 0; path_idx < all_paths.size(); path_idx++) {
                                    auto node = all_paths[path_idx][all_paths[path_idx].size()-1];
                                    if (node.first == previous_contig) {
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
                                            coverages[path[i].first] += min(1.0, (line.size() - pos_begin) / (double)length_of_contigs[path[i].first]);
                                            omp_unset_lock(&coverage_lock);
                                        }
                                    }
                                }
                                
                                if (valid_path_count != 1 && current_path.contigs.size() > 0) {
                                    paths.push_back(current_path);
                                    current_path.contigs.clear();
                                    current_path.orientations.clear();
                                    current_path.start_position_on_contig = length_of_contigs[contig] - pos_in_contig - km;
                                    current_path.end_position_on_contig = length_of_contigs[contig] - pos_in_contig;
                                }
                            }

                            current_path.contigs.push_back(contig);
                            current_path.orientations.push_back(false);
                            if (current_path.start_position_on_contig == -1) {
                                current_path.start_position_on_contig = length_of_contigs[contig] - pos_in_contig - km;
                            }
                            current_path.end_position_on_contig = length_of_contigs[contig] - pos_in_contig;
                            
                            omp_set_lock(&coverage_lock);
                            coverages[contig]+= min(1.0, (line.size()- pos_begin) / (double)length_of_contigs[contig]);
                            omp_unset_lock(&coverage_lock);
                        }
                        
                        int length_of_contig_left = pos_in_contig;
                        if (length_of_contig_left > 10){
                            pos_to_look_at += min((int) (length_of_contig_left*0.8), (int)(line.size() - pos_begin - km - 5));
                        }
                        previous_match = pos_begin;
                        previous_contig = contig;
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
                        thread_output << (p.orientations[i] ? ">" : "<") << p.contigs[i];
                    }
                    thread_output << "\t";
                    
                    int path_length = 0;
                    for (const auto& contig : p.contigs){
                        path_length += length_of_contigs[contig];
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

    add_coverages_to_graph(unitig_graph, coverages);

    cout << "    -> Number of cleanly corrected/aligned reads: " << nb_reads_single_path << " out of " << nb_reads << endl;
}


/**
 * @brief Given a starting position on a contig and an orientation, follow the graph and list all possible contigs and paths
 * 
 * @param linked adjacency list representing the graph structure
 * @param start_contig_name name of the starting contig
 * @param start_orientation true if we start from the right end ('+'), false if from the left end ('-')
 * @param max_length maximum total length of contigs to explore in the paths
 * @param length_of_contigs map from contig name to its length
 * @return vector<vector<pair<string, bool>>> list of paths, each path is a list of (contig_name, orientation)
 */
vector<vector<pair<string, bool>>> list_all_paths_from_contig(unordered_map<string, pair<vector<pair<string, char>>, vector<pair<string,char>>>>& linked, const string& start_contig_name, bool start_orientation, int max_length, unordered_flat_map<string, int>& length_of_contigs, int km){
    
    vector<vector<pair<string, bool>>> all_results;
    all_results.reserve(50); // Reserve space for expected maximum paths
    
    // Recursive exploration function using push/pop
    std::function<void(const string&, char, int, vector<pair<string, bool>>&)> explore;
    explore = [&](const string& current_contig, char current_end, int length_left, vector<pair<string, bool>>& current_path) {
        
        if (length_left <= 0){
            all_results.push_back(current_path);
            return;
        }
        
        // Early termination if too many paths
        if (all_results.size() >= 50){
            return;
        }
        
        // Get neighbors from the appropriate end
        const vector<pair<string, char>>& neighbors = (current_end == 1) ? linked[current_contig].second : linked[current_contig].first;
        
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
            
            // Continue exploration from the opposite end of the neighbor
            char next_end = 1 - neighbor.second;
            int neighbor_length = length_of_contigs[neighbor.first] - km + 1;
            
            explore(neighbor.first, next_end, length_left - neighbor_length, current_path);
            
            // Pop from path (backtrack)
            current_path.pop_back();
            
            // Early termination check
            if (all_results.size() >= 50){
                return;
            }
        }
    };
    
    // Start exploration
    vector<pair<string, bool>> initial_path;
    initial_path.reserve(20); // Reserve space for typical path length
    initial_path.push_back({start_contig_name, start_orientation});
    
    char start_end = start_orientation ? 1 : 0; // if orientation is true ('+'), we start from end 1 (right)
    
    explore(start_contig_name, start_end, max_length, initial_path);
    
    // Return empty if we hit the limit
    if (all_results.size() >= 50){
        return {};
    }
    
    return all_results;
}



void merge_adjacent_contigs_BCALM(std::string gfa_in, std::string gfa_out, int k, std::string path_to_bcalm, std::string path_convertToGFA, std::string path_tmp_folder){
        
        //convert gfa_in to fasta
        string tmp_fasta = path_tmp_folder + "tmp_324.fasta";
        gfa_to_fasta(gfa_in, tmp_fasta);

        //to merge, simply make a unitig graph from bcalm.unitigs.shaved.gfa and then convert it to gfa
        // cout << "Creating shaved unitig graph\n";
        string command_unitig_graph = path_to_bcalm + " -in " + tmp_fasta + " -kmer-size "+std::to_string(k)+" -abundance-min 1 -out "+path_tmp_folder+"tmp_324 > "+path_tmp_folder+"bcalm.log 2>&1";
        auto unitig_graph_ok = system(command_unitig_graph.c_str());
        // cout << "launching unitig graph\n" << command_unitig_graph << endl;
        if (unitig_graph_ok != 0){
            cerr << "ERROR: unitig graph failed in merge_adjacent_contigs_BCALM\n";
            cout << command_unitig_graph << endl;
            exit(1);
        }

        //convert to gfa
        // cout << "Launching convertToGFA\n";
        string convert_command2 = path_convertToGFA + " " + path_tmp_folder+ "tmp_324.unitigs.fa " + gfa_out + " " + std::to_string(k) + " > "+path_tmp_folder+"convertToGFA.log 2>&1";
        auto res = system(convert_command2.c_str());

        //remove tmp files
        // string remove_tmp_files = "rm "+path_tmp_folder+"tmp_324*";
        // system(remove_tmp_files.c_str());
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

    //move the sorted file to the original file
    std::string command = "mv " + gfa + ".tmp " + gfa;
    auto res = system(command.c_str());

}

//recursive function that returns true if it can find a path on one side of the contig of the contig of size 4*k without going through a contig with a coverage more than 20x the coverage of the contig
//three possible outcomes: 0 = nothing overcovered but dead end, 1 = overcovered, 2 = not overcovered path found
/**
 * @brief 
 * 
 * @param contig current contig
 * @param endOfContig endOfContig we arrive from
 * @param linked 
 * @param coverage 
 * @param length_of_contigs 
 * @param k 
 * @param length_left 
 * @param original_coverage 
 * @param not_overcovered set to true if the contig is not overcovered left or right
 * @return int 2: not overcovered path found, 1: overcovered, 0: nothing overcovered but dead end
 */
int explore_neighborhood(string contig, int endOfContig, unordered_map<string, pair<vector<pair<string, char>>, vector<pair<string,char>>>> &linked, unordered_map<string, float>& coverage, unordered_map<string, int>& length_of_contigs, int k, int length_left, double original_coverage, double big_coverage, unordered_map<string, pair<bool, bool>> &not_overcovered){
    
    if (length_left <= 0){

        if (coverage[contig] > big_coverage){
            return 1;
        }
        else if (endOfContig == 1 && linked[contig].first.size() == 0 || endOfContig == 0 && linked[contig].second.size() == 0){
            return 0;
        }
        else{
            return 2;
        }
    }

    bool overcovered_in_neighborhood = false;
    if (endOfContig == 1){ //arrived by the right, go through to the left
        if (linked[contig].first.size() == 0){
            return 0;
        }
        for (auto l: linked[contig].first){
            if (coverage[l.first] > big_coverage){
                overcovered_in_neighborhood = true;
            }
            else if (l.first != contig) { //the condition is so we don't end up in a loop

                //check if we're encountering an already known not overcovered contig
                if (not_overcovered.find(l.first) != not_overcovered.end() && not_overcovered[l.first].first && l.second == 1 && coverage[l.first] < 2*coverage[contig] && coverage[l.first]*2 > original_coverage){
                    return 2;
                }
                else if (not_overcovered.find(l.first) != not_overcovered.end() && not_overcovered[l.first].second && l.second == 0 && coverage[l.first] < 2*coverage[contig] && coverage[l.first]*2 > original_coverage){
                    return 2;
                }

                int res = explore_neighborhood(l.first, l.second, linked, coverage, length_of_contigs, k, length_left - length_of_contigs[l.first] + k - 1, original_coverage, big_coverage, not_overcovered);
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
        if (linked[contig].second.size() == 0){
            return 0;
        }
        for (auto l: linked[contig].second){
            if (coverage[l.first] > big_coverage){
                overcovered_in_neighborhood = true;
            }
            else {

                if (not_overcovered.find(l.first) != not_overcovered.end() && not_overcovered[l.first].first && l.second == 1 && coverage[l.first] < 2*original_coverage && coverage[l.first]*2 > original_coverage){
                    return 2;
                }
                else if (not_overcovered.find(l.first) != not_overcovered.end() && not_overcovered[l.first].second && l.second == 0 && coverage[l.first] < 2*original_coverage && coverage[l.first]*2 > original_coverage){
                    return 2;
                }

                int res = explore_neighborhood(l.first, l.second, linked, coverage, length_of_contigs, k, length_left - length_of_contigs[l.first] + k - 1, original_coverage, big_coverage, not_overcovered);
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

    ifstream input(gfa_in);
    //first go through the gfa and find all the places where an end of contig is connected with two links
    unordered_map<string, pair<vector<pair<string, char>>, vector<pair<string,char>>>> linked;
    unordered_map<string, long int> pos_of_contig_seq_in_file;
    unordered_map<string, float> coverage;
    unordered_map<string, int> length_of_contigs;
    vector<string> list_of_contigs;

    //parse the gfa file
    string line;
    long int pos = 0;
    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            string depth_string = "  ";
            while (depth_string.size()>= 2 && (depth_string.substr(0,2) != "DP" && depth_string.substr(0,2) != "km")){
                string nds;
                ss >> nds;
                depth_string = nds;
            }

            if (depth_string.substr(0,2) != "DP" && depth_string.substr(0,2) != "km"){
                cerr << "ERROR: no depth found for contig " << name << "\n";
                exit(1);
            }
            double depth = std::stof(depth_string.substr(5, depth_string.size()-5));
            depth = std::max(1.0, depth- extra_coverage);
            coverage[name] = depth;
            pos_of_contig_seq_in_file[name] = pos;
            length_of_contigs[name] = sequence.size();
            list_of_contigs.push_back(name);

            if (linked.find(name) == linked.end()){
                linked[name] = {vector<pair<string,char>>(0), vector<pair<string,char>>(0)};
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

            char or1 = (orientation1 == "+" ? 1 : 0);
            char or2 = (orientation2 == "+" ? 0 : 1);
            auto neighbor = make_pair(name2, or2);
            if (orientation1 == "+"){
                if (std::find(linked[name1].second.begin(), linked[name1].second.end(), neighbor) == linked[name1].second.end()){
                    linked[name1].second.push_back(neighbor);
                }
            }
            else{
                if (std::find(linked[name1].first.begin(), linked[name1].first.end(), neighbor) == linked[name1].first.end()){
                    linked[name1].first.push_back(neighbor);
                }
            }

            neighbor = make_pair(name1, or1);
            if (orientation2 == "+"){
                if (std::find(linked[name2].first.begin(), linked[name2].first.end(), neighbor) == linked[name2].first.end()){
                    linked[name2].first.push_back(neighbor);
                }
            }
            else{
                if (std::find(linked[name2].second.begin(), linked[name2].second.end(), neighbor) == linked[name2].second.end()){
                    linked[name2].second.push_back(neighbor);
                }
            }
        }
        pos += line.size() + 1;
    }
    input.close();
    input.open(gfa_in);

    unordered_set<string> to_keep; //kept contigs are the one with a coverage above abundance_min and their necessary neighbors for the contiguity

    //iterative cleaning of the graph

    unordered_map<string, pair<bool, bool>> not_overcovered; //iteratively mark the contigs that are not overcovered left or right
    not_overcovered.reserve(linked.size());
    omp_lock_t lock_contigs_to_keep;
    omp_init_lock(&lock_contigs_to_keep);
    omp_lock_t lock_not_overcovered_contigs;
    omp_init_lock(&lock_not_overcovered_contigs);

    //decide which contig we really want to keep
    #pragma omp parallel for num_threads(num_threads)
    for (int c = 0 ; c < list_of_contigs.size() ; c++){

        string contig = list_of_contigs[c];

        //check if this is a badly covered bubble
        bool bubble = false;
        if (linked[contig].first.size() == 1 && linked[contig].second.size() == 1){

            string neighbor_left = linked[contig].first[0].first;
            char end_of_neighbor_left = linked[contig].first[0].second;
            string neighbor_right = linked[contig].second[0].first;
            char end_of_neighbor_right = linked[contig].second[0].second;

            string other_neighbor_of_contig_left = "";
            if (end_of_neighbor_left == 0 && linked[neighbor_left].first.size() == 2){
                for (auto l: linked[neighbor_left].first){
                    if (l.first != contig){
                        other_neighbor_of_contig_left = l.first;
                    }
                }
            }
            else if (end_of_neighbor_left == 1 && linked[neighbor_left].second.size() == 2){
                for (auto l: linked[neighbor_left].second){
                    if (l.first != contig){
                        other_neighbor_of_contig_left = l.first;
                    }
                }
            }

            string other_neighbor_of_contig_right = "";
            if (end_of_neighbor_right == 0 && linked[neighbor_right].first.size() == 2){
                for (auto l: linked[neighbor_right].first){
                    if (l.first != contig){
                        other_neighbor_of_contig_right = l.first;
                    }
                }
            }
            else if (end_of_neighbor_right == 1 && linked[neighbor_right].second.size() == 2){
                for (auto l: linked[neighbor_right].second){
                    if (l.first != contig){
                        other_neighbor_of_contig_right = l.first;
                    }
                }
            }

            if (other_neighbor_of_contig_left == other_neighbor_of_contig_right && other_neighbor_of_contig_left != "" && 5*coverage[contig] < coverage[other_neighbor_of_contig_left]){
                bubble = true;
            }
            if (single_genome && other_neighbor_of_contig_left == other_neighbor_of_contig_right && other_neighbor_of_contig_left != "" && 2*coverage[contig] < coverage[other_neighbor_of_contig_left]){
                bubble = true;
            }
    
        }

        if (bubble && ((coverage[contig] < abundance_min || abundance_min == -1) && length_of_contigs[contig] < min_length)) //pop it if 1) this is a bubble and 2) coverage less than 5x and the contig is shorter than min_length
        {
            //do nothing, and most importantly, do not add the contig to the to_keep set
        }
        else if (coverage[contig] > abundance_min || length_of_contigs[contig] > min_length){ 
            
            int size_of_neighborhood = 7*k;
            // cout << "launching..\n";
            int overcovered_right = 2;
            if (not_overcovered.find(contig) == not_overcovered.end() || !not_overcovered[contig].second){
                double big_coverage = 20*coverage[contig];
                if (coverage[contig] < 5){ //not very solid, don't make such a fuss about deleting it
                    big_coverage = 3*coverage[contig];
                }
                overcovered_right = explore_neighborhood(contig, 0, linked, coverage, length_of_contigs, k, size_of_neighborhood, coverage[contig], big_coverage, not_overcovered);
            }
            if (overcovered_right == 2){
                omp_set_lock(&lock_not_overcovered_contigs);
                
                if (not_overcovered.find(contig) != not_overcovered.end()){
                    not_overcovered[contig].first = false;
                    not_overcovered[contig].second = false;
                }
                not_overcovered[contig].second = true;

                omp_unset_lock(&lock_not_overcovered_contigs);
            }
            int overcovered_left = 2;
            if (not_overcovered.find(contig) == not_overcovered.end() || !not_overcovered[contig].first){
                double big_coverage = 20*coverage[contig];
                if (coverage[contig] < 5){ //not very solid, don't make such a fuss about deleting it
                    big_coverage = 3*coverage[contig];
                }
                overcovered_left = explore_neighborhood(contig, 1, linked, coverage, length_of_contigs, k, size_of_neighborhood, coverage[contig], big_coverage, not_overcovered);
            }
            if (overcovered_left == 2){
                omp_set_lock(&lock_not_overcovered_contigs);
                if (not_overcovered.find(contig) != not_overcovered.end()){
                    not_overcovered[contig].first = false;
                    not_overcovered[contig].second = false;
                }
                not_overcovered[contig].first = true;
                omp_unset_lock(&lock_not_overcovered_contigs);
            }

            //now decide if the contig is to be kept
            if (overcovered_left == 1 && overcovered_right == 1){ //overcovered on both sides, pop it

            }
            else if (overcovered_left == 0 && overcovered_right == 0 && length_of_contigs[contig] > min_length){ //if the contig has two dead ends, keep it under conditiosn that it is long enough
                omp_set_lock(&lock_contigs_to_keep);
                to_keep.insert(contig);
                omp_unset_lock(&lock_contigs_to_keep);
            }
            else if (overcovered_left == 2 || overcovered_right == 2){ //if the contig is not overcovered on one side, keep it (and it also passed the abundance_min or the min_length threshold)
                omp_set_lock(&lock_contigs_to_keep);
                to_keep.insert(contig);
                omp_unset_lock(&lock_contigs_to_keep);
            }
            else if (overcovered_left == 1 && overcovered_right == 0 || overcovered_left == 0 && overcovered_right == 1){ //this means that this is a tip
                //do nothing, and most importantly, do not add the contig to the to_keep set
            }

            //if single genome mode is on, delete badly covered dead ends
            if (single_genome) {
                if (linked[contig].first.size() == 0 && linked[contig].second.size() > 0  ){
                    float max_neighbor_coverage = 0;
                    for (auto neighbor : linked[contig].second) {
                        max_neighbor_coverage = std::max(max_neighbor_coverage, coverage[neighbor.first]);
                    }
                    if (coverage[contig] * 2 < max_neighbor_coverage && length_of_contigs[contig] < 2*min_length) {
                        omp_set_lock(&lock_contigs_to_keep);
                        to_keep.erase(contig);
                        omp_unset_lock(&lock_contigs_to_keep);
                    }
                }
                if (linked[contig].first.size() > 0 && linked[contig].second.size() == 0  ){
                    float max_neighbor_coverage = 0;
                    for (auto neighbor : linked[contig].first) {
                        max_neighbor_coverage = std::max(max_neighbor_coverage, coverage[neighbor.first]);
                    }
                    if (coverage[contig] * 2 < max_neighbor_coverage && length_of_contigs[contig] < 2*min_length) {
                        omp_set_lock(&lock_contigs_to_keep);
                        to_keep.erase(contig);
                        omp_unset_lock(&lock_contigs_to_keep);
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

            string contig = list_of_contigs[c];

            if (to_keep.find(contig) != to_keep.end()){
                //make sure the contig has at least one neighbor left and right (if not, take the one with the highest coverage) 
                //this is the only way to keep contigs that are below abundance_min
                float best_coverage = 0;
                string best_contig = "";
                bool at_least_one_neighbor = false;
                for (auto l: linked[contig].second){
                    if (coverage[l.first] > best_coverage){
                        best_coverage = coverage[l.first];
                        best_contig = l.first;
                    }
                    if (to_keep.find(l.first) != to_keep.end()){
                        at_least_one_neighbor = true;
                    }
                }
                if (best_contig != "" && !at_least_one_neighbor){
                    if (to_keep.find(best_contig) == to_keep.end()){
                        omp_set_lock(&lock_contigs_to_keep);
                        to_keep.insert(best_contig);
                        omp_unset_lock(&lock_contigs_to_keep);
                        number_of_edits++;
                    }
                }

                best_coverage = 0;
                best_contig = "";
                at_least_one_neighbor = false;
                for (auto l: linked[contig].first){
                    if (coverage[l.first] > best_coverage){
                        best_coverage = coverage[l.first];
                        best_contig = l.first;
                    }
                    if (to_keep.find(l.first) != to_keep.end()){
                        at_least_one_neighbor = true;
                    }
                }
                if (best_contig != "" && !at_least_one_neighbor){
                    if (to_keep.find(best_contig) == to_keep.end()){
                        omp_set_lock(&lock_contigs_to_keep);
                        to_keep.insert(best_contig);
                        omp_unset_lock(&lock_contigs_to_keep);
                        number_of_edits++;
                    }
                }
            }
        }
    }
    omp_destroy_lock(&lock_contigs_to_keep);
    omp_destroy_lock(&lock_not_overcovered_contigs);

    //now wirte the gfa file without the contigs to remove
    input.close();
    input.open(gfa_in);
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
            if (to_keep.find(name) != to_keep.end()){
                out << "S\t" << name << "\t" << sequence << "\tLN:i:" << sequence.size() << "\tkm:f:" << coverage[name] << "\n";
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

            if (to_keep.find(name1) != to_keep.end() && to_keep.find(name2) != to_keep.end()){
                out << line << "\n";
            }
        }
    }
    out.close();
}

/**
 * @brief Function that takes as input a graph, cleans it for building long contigs to ouptut for next k. The goal is to improve contiguity at next step
 * 
 * @param gfa_in 
 * @param gfa_out 
 * @param k
 * @param extra_coverage //to retreat to the coverage because it comes from extra contigs added to the reads from previous assembly rounds
 * @param num_threads
 */
void trim_graph_for_next_k(string gfa_in, string gfa_out,  int k, int extra_coverage, int num_threads){

    ifstream input(gfa_in);
    //first go through the gfa and find all the places where an end of contig is connected with two links
    unordered_map<string, pair<vector<pair<string, char>>, vector<pair<string,char>>>> linked;
    unordered_map<string, pair<int,int>> number_of_links; //number of links on each side of the contig

    unordered_map<string, long int> pos_of_contig_seq_in_file;
    unordered_map<string, float> coverage;
    unordered_map<string, int> length_of_contigs;
    vector<string> list_of_contigs;

    //parse the gfa file
    string line;
    long int pos = 0;
    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            string depth_string = "  ";
            while (depth_string.size()>= 2 && (depth_string.substr(0,2) != "DP" && depth_string.substr(0,2) != "km")){
                string nds;
                ss >> nds;
                depth_string = nds;
            }

            if (depth_string.substr(0,2) != "DP" && depth_string.substr(0,2) != "km"){
                cerr << "ERROR: no depth found for contig " << name << "\n";
                exit(1);
            }
            double depth = std::stof(depth_string.substr(5, depth_string.size()-5));
            depth = std::max(1.0, depth- extra_coverage);
            coverage[name] = depth;
            pos_of_contig_seq_in_file[name] = pos;
            length_of_contigs[name] = sequence.size();
            list_of_contigs.push_back(name);

            if (linked.find(name) == linked.end()){
                linked[name] = {vector<pair<string,char>>(0), vector<pair<string,char>>(0)};
                number_of_links[name] = {0, 0};
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

            char or1 = (orientation1 == "+" ? 1 : 0);
            char or2 = (orientation2 == "+" ? 0 : 1);
            auto neighbor = make_pair(name2, or2);
            if (orientation1 == "+"){
                if (std::find(linked[name1].second.begin(), linked[name1].second.end(), neighbor) == linked[name1].second.end()){
                    linked[name1].second.push_back(neighbor);
                    number_of_links[name1].second++;
                }
            }
            else{
                if (std::find(linked[name1].first.begin(), linked[name1].first.end(), neighbor) == linked[name1].first.end()){
                    linked[name1].first.push_back(neighbor);
                    number_of_links[name1].first++;
                }
            }

            neighbor = make_pair(name1, or1);
            if (orientation2 == "+"){
                if (std::find(linked[name2].first.begin(), linked[name2].first.end(), neighbor) == linked[name2].first.end()){
                    linked[name2].first.push_back(neighbor);
                    number_of_links[name2].first++;
                }
            }
            else{
                if (std::find(linked[name2].second.begin(), linked[name2].second.end(), neighbor) == linked[name2].second.end()){
                    linked[name2].second.push_back(neighbor);
                    number_of_links[name2].second++;
                }
            }
        }
        pos += line.size() + 1;
    }
    input.close();
    input.open(gfa_in);


    std::set<pair<pair<string,char>,pair<string,char>>> links_to_erase; //links to erase

    //cleaning of the graph

    omp_lock_t lock_to_erase;
    omp_init_lock(&lock_to_erase);

    int nb_changes = 11;
    while (nb_changes > 10){
        nb_changes = 0;

        //find the tips in parallel
        #pragma omp parallel for num_threads(num_threads)
        for (int c = 0 ; c < list_of_contigs.size() ; c++){

            string contig = list_of_contigs[c];


            if (number_of_links[contig].first*number_of_links[contig].second == 0 && number_of_links[contig].first + number_of_links[contig].second == 1){

                char end_of_contig = (number_of_links[contig].first == 1 ? 0 : 1);
                string neighbor = (end_of_contig == 0 ? linked[contig].first[0].first : linked[contig].second[0].first);
                char end_of_neighbor = (end_of_contig == 0 ? linked[contig].first[0].second : linked[contig].second[0].second);
                auto& neighbors_of_neighbor = (end_of_neighbor == 0 ? linked[neighbor].first : linked[neighbor].second);

                for (auto l: neighbors_of_neighbor){
                    if (l.first == contig){
                        continue;
                    }
                    char end_of_neighbor_neighbor = l.second;
                    int number_of_links_other_side_neighbor_of_neighbor = (end_of_neighbor_neighbor == 0 ? linked[l.first].second.size() : linked[l.first].first.size());

                    if (coverage[l.first] > 2*coverage[contig] || length_of_contigs[contig] < k + 10 || number_of_links_other_side_neighbor_of_neighbor > 0){
                        omp_set_lock(&lock_to_erase);
                        links_to_erase.insert(make_pair(make_pair(contig, end_of_contig), make_pair(neighbor, end_of_neighbor)));
                        links_to_erase.insert(make_pair(make_pair(neighbor, end_of_neighbor), make_pair(contig, end_of_contig)));
                        if (end_of_contig == 0){
                            number_of_links[contig].first--;
                        }
                        else{
                            number_of_links[contig].second--;
                        }
                        if (end_of_neighbor == 0){
                            number_of_links[neighbor].first--;
                        }
                        else{
                            number_of_links[neighbor].second--;
                        }
                        nb_changes++;

                        omp_unset_lock(&lock_to_erase);
                    }
                }
            }
        }
    }

    omp_destroy_lock(&lock_to_erase);

    //now wirte the gfa file without the links to remove
    input.close();
    input.open(gfa_in);
    ofstream out(gfa_out);
    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            out << line << "\n";
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

            char or1 = (orientation1 == "+" ? 1 : 0);
            char or2 = (orientation2 == "+" ? 0 : 1);
            if (links_to_erase.find(make_pair(make_pair(name1, or1), make_pair(name2, or2))) == links_to_erase.end()){
                out << line << "\n";
            }
        }
    }
    out.close();
}


/**
 * @brief In small bubbles, take only one side and discard the other
 * 
 * @param gfa_in 
 * @param length_of_longest_read only pop bubbles that are between two contigs that are at least this length and at most this length/2 apart
 * @param gfa_out 
 */
void cut_links_for_contiguity(std::string gfa_in, std::string gfa_out){

    //load the graph
    ifstream input(gfa_in);
    unordered_map<string, pair<vector<pair<string, char>>, vector<pair<string,char>>>> linked;
    unordered_map<string, long int> pos_of_contig_seq_in_file;
    unordered_map<string, float> coverage;
    unordered_map<string, int> length_of_contigs;

    string line;
    long int pos = 0;
    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            string depth_string = "  ";
            while (depth_string.size()>= 2 && (depth_string.substr(0,2) != "DP" && depth_string.substr(0,2) != "km")){
                string nds;
                ss >> nds;
                depth_string = nds;
            }

            if (depth_string.substr(0,2) != "DP" && depth_string.substr(0,2) != "km"){
                cerr << "ERROR: no depth found for contig " << name << "\n";
                exit(1);
            }
            float depth = std::stof(depth_string.substr(5, depth_string.size()-5));
            coverage[name] = depth;
            pos_of_contig_seq_in_file[name] = pos;
            length_of_contigs[name] = sequence.size();

            if (linked.find(name) == linked.end()){
                linked[name] = {vector<pair<string,char>>(0), vector<pair<string,char>>(0)};
            }
        }
    }
    input.close();
    input.open(gfa_in);
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

            if (linked.find(name1) == linked.end() || linked.find(name2) == linked.end()){
                cerr << "ERROR: contig not found in linked in pop_bubbles: ";
                cerr << name1 << " " << name2 << "\n";
                exit(1);
            }

            char or1 = (orientation1 == "+" ? 1 : 0);
            char or2 = (orientation2 == "+" ? 0 : 1);
            auto neighbor = make_pair(name2, or2);
            if (orientation1 == "+"){
                if (std::find(linked[name1].second.begin(), linked[name1].second.end(), neighbor) == linked[name1].second.end()){
                    linked[name1].second.push_back(neighbor);
                }
            }
            else{
                if (std::find(linked[name1].first.begin(), linked[name1].first.end(), neighbor) == linked[name1].first.end()){
                    linked[name1].first.push_back(neighbor);
                }
            }

            neighbor = make_pair(name1, or1);
            if (orientation2 == "+"){
                if (std::find(linked[name2].first.begin(), linked[name2].first.end(), neighbor) == linked[name2].first.end()){
                    linked[name2].first.push_back(neighbor);
                }
            }
            else{
                if (std::find(linked[name2].second.begin(), linked[name2].second.end(), neighbor) == linked[name2].second.end()){
                    linked[name2].second.push_back(neighbor);
                }
            }
        }
        pos += line.size() + 1;
    }
    input.close();

    //now cut the links that are not good for contiguity

    std::set<pair<pair<string,char>,pair<string,char>>> links_to_delete;
    for (auto c: linked){
        string contig_name = c.first;
        for (char end = 0 ; end < 2 ; end++){
            vector<pair<string, char>>& neighbors = linked[contig_name].first;
            if (end == 1){
                neighbors = linked[contig_name].second;
            }
            if (neighbors.size() > 1){
                float max_coverage = 0;
                bool all_neighbors_are_dead_ends = true;
                for (auto& neighbor : neighbors){
                    if (coverage[neighbor.first] > max_coverage){
                        max_coverage = coverage[neighbor.first];
                    }
                    if (linked[neighbor.first].first.size() > 0 && linked[neighbor.first].second.size() > 0){
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
                    //if the neighbor is a dead end cut it
                    if (!all_neighbors_are_dead_ends &&
                            (linked[neighbor.first].first.size() == 0 || linked[neighbor.first].second.size() == 0)){
                        links_to_delete.insert({neighbor, {contig_name,end}});
                        links_to_delete.insert({{contig_name,end}, neighbor});
                    }
                }
            }
        }
        
    }

    //now write the gfa file without the links to delete
    input.open(gfa_in);
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

            if (links_to_delete.find({{name1, (orientation1 == "+" ? 1 : 0)}, {name2, (orientation2 == "+" ? 0 : 1)}}) == links_to_delete.end()){
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
void trim_tips_isolated_contigs_and_bubbles(std::string gfa_in, int min_coverage, int min_length, std::string gfa_out){
    //load the graph
    ifstream input(gfa_in);
    unordered_map<string, pair<vector<pair<string, char>>, vector<pair<string,char>>>> linked;
    unordered_map<string, long int> pos_of_contig_seq_in_file;
    unordered_map<string, float> coverage;
    unordered_map<string, int> length_of_contigs;

    string line;
    long int pos = 0;
    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            string depth_string = "  ";
            while (depth_string.size()>= 2 && (depth_string.substr(0,2) != "DP" && depth_string.substr(0,2) != "km")){
                string nds;
                ss >> nds;
                depth_string = nds;
            }

            if (depth_string.substr(0,2) != "DP" && depth_string.substr(0,2) != "km"){
                cerr << "ERROR: no depth found for contig " << name << "\n";
                exit(1);
            }
            float depth = std::stof(depth_string.substr(5, depth_string.size()-5));
            coverage[name] = depth;
            pos_of_contig_seq_in_file[name] = pos;
            length_of_contigs[name] = sequence.size();

            if (linked.find(name) == linked.end()){
                linked[name] = {vector<pair<string,char>>(0), vector<pair<string,char>>(0)};
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

            if (linked.find(name1) == linked.end() || linked.find(name2) == linked.end()){
                cerr << "ERROR: contig not found in linked in pop_bubbles: ";
                cerr << name1 << " " << name2 << "\n";
                exit(1);
            }

            char or1 = (orientation1 == "+" ? 1 : 0);
            char or2 = (orientation2 == "+" ? 0 : 1);
            auto neighbor = make_pair(name2, or2);
            if (orientation1 == "+"){
                if (std::find(linked[name1].second.begin(), linked[name1].second.end(), neighbor) == linked[name1].second.end()){
                    linked[name1].second.push_back(neighbor);
                }
            }
            else{
                if (std::find(linked[name1].first.begin(), linked[name1].first.end(), neighbor) == linked[name1].first.end()){
                    linked[name1].first.push_back(neighbor);
                }
            }

            neighbor = make_pair(name1, or1);
            if (orientation2 == "+"){
                if (std::find(linked[name2].first.begin(), linked[name2].first.end(), neighbor) == linked[name2].first.end()){
                    linked[name2].first.push_back(neighbor);
                }
            }
            else{
                if (std::find(linked[name2].second.begin(), linked[name2].second.end(), neighbor) == linked[name2].second.end()){
                    linked[name2].second.push_back(neighbor);
                }
            }
        }
        pos += line.size() + 1;
    }
    input.close();
    
    //now trim the tips, isolated contigs and bubbles with a coverage below min_coverage and a length below min_length (bubbles need to have coverage 1)
    std::set<string> contigs_to_remove;
    std::set<std::pair<string, string>> links_to_detach; //contigs that look valid but probably reduce contiguity, detach them from the graph
    for (auto c: linked){
        string contig_name = c.first;
        //remove tip
        if (linked[contig_name].first.size() == 0 || linked[contig_name].second.size() == 0){
            if (coverage[contig_name] < min_coverage && (length_of_contigs[contig_name] < min_length || coverage[contig_name]==1)){
                contigs_to_remove.insert(contig_name);
            }
            else { //the contig is valid but detaching it may improve the contiguity
                if (linked[contig_name].first.size() > 0) {
                    for (auto l : linked[contig_name].first) {
                        if (coverage[l.first] > 2 * coverage[contig_name]) {
                            links_to_detach.insert({contig_name, l.first});
                        }
                    }
                }
                if (linked[contig_name].second.size() > 0) {
                    for (auto l : linked[contig_name].second) {
                        if (coverage[l.first] > 2 * coverage[contig_name]) {
                            links_to_detach.insert({contig_name, l.first});
                        }
                    }
                }
            }
        }
        //remove bubble
        if (linked[contig_name].first.size() == 1 && linked[contig_name].second.size() == 1 && coverage[contig_name] == 1) {
            string neighbor_left = linked[contig_name].first[0].first;
            char end_of_neighbor_left = linked[contig_name].first[0].second;
            string neighbor_right = linked[contig_name].second[0].first;
            char end_of_neighbor_right = linked[contig_name].second[0].second;

            string other_neighbor_of_contig_left = "";
            if (end_of_neighbor_left == 0 && linked[neighbor_left].first.size() == 2) {
                for (auto l : linked[neighbor_left].first) {
                    if (l.first != contig_name) {
                        other_neighbor_of_contig_left = l.first;
                    }
                }
            } else if (end_of_neighbor_left == 1 && linked[neighbor_left].second.size() == 2) {
                for (auto l : linked[neighbor_left].second) {
                    if (l.first != contig_name) {
                        other_neighbor_of_contig_left = l.first;
                    }
                }
            }

            string other_neighbor_of_contig_right = "";
            if (end_of_neighbor_right == 0 && linked[neighbor_right].first.size() == 2) {
                for (auto l : linked[neighbor_right].first) {
                    if (l.first != contig_name) {
                        other_neighbor_of_contig_right = l.first;
                    }
                }
            } else if (end_of_neighbor_right == 1 && linked[neighbor_right].second.size() == 2) {
                for (auto l : linked[neighbor_right].second) {
                    if (l.first != contig_name) {
                        other_neighbor_of_contig_right = l.first;
                    }
                }
            }

            if (other_neighbor_of_contig_left == other_neighbor_of_contig_right && other_neighbor_of_contig_left != "") {
                contigs_to_remove.insert(contig_name);
            }
        }
        //if contig attached among contigs of higher coverage, detach
        if (linked[contig_name].first.size() > 0 && linked[contig_name].second.size() > 0) {
            // Check if all neighbors have significantly higher coverage
            bool all_neighbors_higher_coverage = true;
            for (auto l : linked[contig_name].first) {
                if (coverage[l.first] < 2 * coverage[contig_name]) {
                    all_neighbors_higher_coverage = false;
                    break;
                }
            }
            if (all_neighbors_higher_coverage) {
                for (auto l : linked[contig_name].second) {
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
                for (auto l : linked[contig_name].first) {
                    string neighbor = l.first;
                    char neighbor_end = l.second;
                    auto& neighbor_links = (neighbor_end == 0) ? linked[neighbor].first : linked[neighbor].second;
                    
                    bool has_better_alternative = false;
                    for (auto alt_link : neighbor_links) {
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
                    for (auto l : linked[contig_name].second) {
                        string neighbor = l.first;
                        char neighbor_end = l.second;
                        auto& neighbor_links = (neighbor_end == 0) ? linked[neighbor].first : linked[neighbor].second;
                        
                        bool has_better_alternative = false;
                        for (auto alt_link : neighbor_links) {
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
                    for (auto l : linked[contig_name].first) {
                        links_to_detach.insert({contig_name, l.first});
                    }
                    for (auto l : linked[contig_name].second) {
                        links_to_detach.insert({contig_name, l.first});
                    }
                }
            }
        }
    }

    //now write the gfa file without the contigs to remove
    input.open(gfa_in);
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
            if (contigs_to_remove.find(name) == contigs_to_remove.end()){
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

            if (contigs_to_remove.find(name1) == contigs_to_remove.end() 
                && contigs_to_remove.find(name2) == contigs_to_remove.end()
                && links_to_detach.find({name1, name2}) == links_to_detach.end() 
                && links_to_detach.find({name2, name1}) == links_to_detach.end()){
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
    while (getline(gfa, line)){
        if (line[0] == 'S'){
            long int pos_in_file = (long int) gfa.tellg() - line.size() - 1;
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
                segments.push_back(s);
            }
            else{
                Segment s(name, segments.size(), vector<pair<vector<pair<int,int>>, vector<string>>>(2), pos_in_file, seq, seq.size(), coverage);
                segment_IDs[name] = s.ID;
                segments.push_back(s);
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

            int ID1 = segment_IDs[name1];
            int ID2 = segment_IDs[name2];

            //check that the link did not already exist
            bool already_exists = false;
            for (pair<int,int> link : segments[ID1].links[end1].first){
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

    set<int> already_looked_at_segments; //old IDs of segments that have already been looked at and merged (don't want to merge them twice)
    unordered_map<pair<int,int>,pair<int,int>> old_ID_to_new_ID; //associates (old_id, old end) with (new_id, new_end)
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
        Segment old_seg = old_segments[seg_idx];

        if (already_looked_at_segments.find(old_seg.ID) != already_looked_at_segments.end()){
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
                int num_matches = std::stoi(cigar.substr(0, cigar.find_first_of("M")));
                all_seqs[thread_num].push_back(seq.substr(num_matches, seq.size()-num_matches));
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
                int num_matches = std::stoi(cigar.substr(0, cigar.find_first_of("M")));
                all_seqs[thread_num].push_back(seq.substr(num_matches, seq.size()-num_matches));
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
            if (already_looked_at_segments.find(old_seg.ID) != already_looked_at_segments.end() || already_looked_at_segments.find(other_end_of_merged_contig_ID) != already_looked_at_segments.end()){
                thread_safe = false;
            }
            else{
                for (int ID : all_IDs[thread_num]){
                    already_looked_at_segments.insert(ID);
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
                for (string name : all_names[thread_num]){
                    new_name += name + "_";
                }
                new_name = new_name.substr(0, new_name.size()-1);
                string new_seq = "";
                for (string seq : all_seqs[thread_num]){
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

                new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), new_length, new_coverage));
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
                for (string name : all_names[thread_num]){
                    new_name += name + "_";
                }
                new_name = new_name.substr(0, new_name.size()-1);
                string new_seq = "";
                for (string seq : all_seqs[thread_num]){
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

                new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), new_length, new_coverage));
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
    for (Segment old_seg : old_segments){
        if (already_looked_at_segments.find(old_seg.ID) == already_looked_at_segments.end()){
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
                already_looked_at_segments.insert(current_ID);
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
                int num_matches = std::stoi(cigar.substr(0, cigar.find_first_of("M")));
                all_seqs.push_back(seq.substr(num_matches, seq.size()-num_matches));
                all_coverages.push_back(old_segments[current_ID].get_coverage());
                all_lengths.push_back(old_segments[current_ID].get_length());
            }
            if (old_segments[current_ID].links[current_end].first[0].first != old_seg.ID){
                circular_as_expected = false;
            }
            if (circular_as_expected){
                already_looked_at_segments.insert(current_ID);
                string new_name = "";
                for (string name : all_names){
                    new_name += name + "_";
                }
                new_name = new_name.substr(0, new_name.size()-1);
                string new_seq = "";
                for (string seq : all_seqs){
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
                new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), new_length, new_coverage));
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

    //now add the links in the new segments
    for (pair<pair<pair<int,int>, pair<int,int>>, string> link : links_to_add){
        new_segments[old_ID_to_new_ID[link.first.first].first].links[old_ID_to_new_ID[link.first.first].second].first.push_back(old_ID_to_new_ID[link.first.second]);
        new_segments[old_ID_to_new_ID[link.first.first].first].links[old_ID_to_new_ID[link.first.first].second].second.push_back(link.second);
    }
}

void output_graph(string gfa_output, string gfa_input, vector<Segment> &segments){
    ofstream gfa(gfa_output);
    for (Segment s : segments){
        if (s.name != "delete_me"){
            gfa << "S\t" << s.name << "\t" << s.get_seq(gfa_input) << "\tDP:f:" << s.get_coverage() <<  "\n";
        }
    }
    int nb_segments_outputted = 0;
    for (Segment s : segments){
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

