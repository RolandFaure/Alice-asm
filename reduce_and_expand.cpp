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
#include <omp.h> //for efficient parallelization
#include <filesystem>
// #include <zlib.h>  // Include the gzstream header
#include <set>


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
 * @brief 
 * 
 * @param reads_file 
 * @param assemblyFile 
 * @param context_length 
 * @param compression 
 * @param km 
 * @param kmers associates a compressed kmer with its positions (central and full) in the kmers_file file
 * @param kmers_file files with the kmers and their reverse complement
 * @param num_threads 
 * @param homopolymer_compression 
 */
void go_through_the_reads_again_and_index_interesting_kmers(string reads_file, 
    string assemblyFile, 
    int order, 
    int compression, 
    int km, 
    std::unordered_set<uint64_t> &central_kmers_in_assembly,
    std::unordered_set<uint64_t> &full_kmers_in_assembly,
    unordered_map<uint64_t, pair<unsigned long long,unsigned long long>>& kmers,
    string central_kmers_file, 
    string full_kmers_file,
    int num_threads, 
    bool homopolymer_compression){

    //empty kmer_files
    std::ofstream out(central_kmers_file);
    out.close();
    out.open(full_kmers_file);
    out.close();

    unordered_map<uint64_t, unordered_map<string, int>> kmer_count; //associates a kmer to a map of all potential uncompressed kmers it corresponds to and their count
    unordered_set<uint64_t> confirmed_kmers;

    ifstream input2(reads_file, std::ios::binary | std::ios::ate);
    std::streamoff file_size = input2.tellg();
    unsigned long long int size_of_chunk = 100000000;
    input2.close();

    int num_full_kmers = 0;

    //go through the fasta file and compute the hash of all the kmer using homecoded ntHash
    int seq_num = 0;
    long long output_limit = 0;
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int chunk = 0 ; chunk <= file_size/size_of_chunk ; chunk++){

        unordered_map<uint64_t, pair<pair<string,string>, bool>> kmers_to_output; //kmers to write in the kmers file from this chunk: associates a kmer to central seq, full seq and and whether it is a confirmed hit when this record was written

        std::ifstream input(reads_file);
        input.seekg(chunk*size_of_chunk);
    
        std::string line;
        bool next_line_is_seq = false;
        string read_name;

        while (std::getline(input, line)){

            if (line[0] == '>')
            {
                if (seq_num > output_limit && omp_get_thread_num() == 0){
                    #pragma omp critical
                    {
                        //display the date and time
                        time_t now = time(0);
                        tm *ltm = localtime(&now);
                        cout << "[" << 1 + ltm->tm_mday << "/" << 1 + ltm->tm_mon << "/" << 1900 + ltm->tm_year << " " << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec << "]" << " Processed " << seq_num << " reads" << endl;
                        output_limit += 50000;
                    }
                    // cout << "nanaaaammmmma " << line << endl;
                }
                seq_num++;
                next_line_is_seq = true;
                read_name = line;
            }
            else if (next_line_is_seq){
                //let's launch the foward and reverse rolling hash
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
                uint64_t hash_foward_compressed, hash_reverse_compressed; //rolling hash of the compressed kmer
                string compressed_read = "";
                compressed_read.reserve(line.size()/compression);

                while (roll(hash_foward, hash_reverse, order, line, pos_end, pos_begin, pos_middle, homopolymer_compression)){
                    if (number_of_hashed_bases >= order){

                        if (line[pos_end] != 'A' && line[pos_end] != 'C' && line[pos_end] != 'G' && line[pos_end] != 'T'){
                            //then finish outputting the line and create a new one
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

                            //roll the compressed kmer
                            roll(hash_foward_compressed, hash_reverse_compressed, km, compressed_read, pos_end_compressed, pos_begin_compressed, pos_middle_compressed, false);
                            
                            positions_sampled.push_back(pos_middle);

                            if (positions_sampled.size() >= km){

                                bool foward_central_kmer_in_assembly = central_kmers_in_assembly.find(hash_foward_compressed) != central_kmers_in_assembly.end();
                                bool reverse_central_kmer_in_assembly = central_kmers_in_assembly.find(hash_reverse_compressed) != central_kmers_in_assembly.end();
                                bool central_kmer_in_assembly = foward_central_kmer_in_assembly || reverse_central_kmer_in_assembly;
                                bool foward_full_kmer_in_assembly = full_kmers_in_assembly.find(hash_foward_compressed) != full_kmers_in_assembly.end();
                                bool reverse_full_kmer_in_assembly = full_kmers_in_assembly.find(hash_reverse_compressed) != full_kmers_in_assembly.end();
                                bool full_kmer_in_assembly = foward_full_kmer_in_assembly || reverse_full_kmer_in_assembly;

                                if (central_kmer_in_assembly || full_kmer_in_assembly){
                                
                                    uint64_t canonical_hash = std::min(hash_foward_compressed, hash_reverse_compressed);
                                    
                                    if (confirmed_kmers.find(canonical_hash) == confirmed_kmers.end()){
                                    
                                        string canonical_seq = "", central_seq = "", reverse_central_seq = "";
                                        if (central_kmer_in_assembly){
                                            central_seq = line.substr(positions_sampled[positions_sampled.size()-km+10], positions_sampled[positions_sampled.size()-1-10] - positions_sampled[positions_sampled.size()-km+10]+1);
                                            reverse_central_seq = reverse_complement(central_seq);
                                            canonical_seq = min(central_seq, reverse_central_seq);
                                        }
                                            
                                        string full_seq="", reverse_full_seq="";
                                        if (full_kmer_in_assembly){
                                            auto begin = std::max(0,positions_sampled[positions_sampled.size()-km] - order);
                                            auto end = std::min((int)line.size()-1, positions_sampled[positions_sampled.size()-1] + order);
                                            full_seq = line.substr(begin, end - begin +1);
                                            reverse_full_seq = reverse_complement(full_seq);
                                            canonical_seq = min(full_seq, reverse_full_seq);
                                        }

                                        #pragma omp critical
                                        {

                                            if (kmer_count.find(canonical_hash) == kmer_count.end()){
                                                kmer_count[canonical_hash] = {};

                                                if (foward_central_kmer_in_assembly || foward_full_kmer_in_assembly){
                                                    kmers_to_output[hash_foward_compressed] = {{central_seq, full_seq}, false};
                                                }
                                                if (reverse_central_kmer_in_assembly || reverse_full_kmer_in_assembly){
                                                    kmers_to_output[hash_reverse_compressed] = {{reverse_central_seq, reverse_full_seq}, false};
                                                }
                                            }

                                            if (kmer_count[canonical_hash].find(canonical_seq) == kmer_count[canonical_hash].end()){
                                                kmer_count[canonical_hash][canonical_seq] = 0;
                                            }
                                            kmer_count[canonical_hash][canonical_seq]++;
                                        
                                            if (kmer_count[canonical_hash][canonical_seq] >  3){
                                                if (foward_central_kmer_in_assembly || foward_full_kmer_in_assembly){
                                                    kmers_to_output[hash_foward_compressed] = {{central_seq, full_seq}, true};
                                                }
                                                if (reverse_central_kmer_in_assembly || reverse_full_kmer_in_assembly){
                                                    kmers_to_output[hash_reverse_compressed] = {{reverse_central_seq, reverse_full_seq}, true};
                                                }
                                                confirmed_kmers.insert(canonical_hash);
                                                kmer_count[canonical_hash].clear();
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    number_of_hashed_bases++;
                }
                //if we went beyond chunk+1 * size_of_chunk, we can stop
                if (input.tellg() > (chunk+1)*size_of_chunk){
                    break;
                }
            } 
        }
        input.close();

        //now write the kmers to the file
        #pragma omp critical
        {
            //open kmer file
            std::ofstream out_central(central_kmers_file, std::ios::app);
            std::ofstream out_full(full_kmers_file, std::ios::app);
            for (auto it = kmers_to_output.begin(); it != kmers_to_output.end(); it++){
                if (it->second.second || confirmed_kmers.find(it->first) == confirmed_kmers.end()){ //else, it means that the kmer was not confirmed at the point where it was introduced but has been confirmed since (probably on another thread), let the confirmed kmer be written by the thread that confirmed it
                    unsigned long long position_central= 1;
                    unsigned long long position_full = 1;
                    if (it->second.first.first != ""){
                        position_central = out_central.tellp();
                        out_central << it->second.first.first << "\n";
                    }
                    if (it->second.first.second != ""){
                        position_full = out_full.tellp();
                        out_full << it->second.first.second << "\n";
                    }
                    kmers[it->first] = {position_central, position_full};
                }
            }
            out_central.close();
            out_full.close();
        }
    }
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
void expand_or_list_kmers_needed_for_expansion(string mode, string asm_reduced, int km, std::unordered_set<uint64_t> &central_kmers_needed, std::unordered_set<uint64_t> &full_kmers_needed, string central_kmers_file, string full_kmers_file, unordered_map<uint64_t, pair<unsigned long long,unsigned long long>>& kmers, string output){
    
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
                    left_seq[contig1] = last_10[contig2];
                    right_seq[contig2] = first_10[contig1];
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
            string expanded_sequence = "";
            int i = 0;
            int length_of_central_kmers = 1 ;
            for (i = 0; i <= sequence.size()-km; i+= km-20-1){ //-20 because we only take the central part of each kmer
                string kmer = sequence.substr(i, km);
                uint64_t hash_foward_kmer = hash_string(km, kmer, false);
                if (mode == "index"){
                    central_kmers_needed.insert(hash_foward_kmer);
                }
                else{
                    if (kmers.find(hash_foward_kmer) != kmers.end() && kmers[hash_foward_kmer].first != 1){
                        //retrieve the central and full sequence from the kmers file
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
                        // cout << "WARNING (code 743) missing kmer " << kmer << "\n";
                        number_of_missing_kmers++;
                        expanded_sequence += string(length_of_central_kmers, 'N');
                    }
                }
            }
            
            // create the beginning of the sequence if there was no left extension
            if (left_seq.find(name) == left_seq.end()){
                string first_kmer = sequence.substr(0, km);
                uint64_t hash_foward_kmer = hash_string(km, first_kmer, false);
                if (mode == "index"){
                    full_kmers_needed.insert(hash_foward_kmer);
                }
                else {
                    if (kmers.find(hash_foward_kmer) != kmers.end() && kmers[hash_foward_kmer].second != 1){

                        //retrieve the full sequence from the kmers file
                        full_kmers_input.seekg(kmers[hash_foward_kmer].second);
                        std::getline(full_kmers_input, line2);
                        string full_kmer = line2;

                        string beginning_of_seq = full_kmer;

                        //compute the overlap
                        int overlap = beginning_of_seq.size();
                        string exp_start = expanded_sequence.substr(0, std::min( (int) expanded_sequence.size(), std::min(30, overlap)));
                        while (overlap > 0 && beginning_of_seq.substr(beginning_of_seq.size()-overlap, exp_start.size()) != exp_start){
                            overlap--;
                            if (overlap < exp_start.size()){
                                exp_start = expanded_sequence.substr(0, overlap);
                            }
                        }
                        expanded_sequence = beginning_of_seq.substr(0, beginning_of_seq.size()-overlap) + expanded_sequence;
                    }
                    else{
                        // cout << "WARNING (code 744) missing kmer " << first_kmer << "\n";
                        number_of_missing_kmers++;  
                    }
                }
            }

            // finish the sequence
            string last_kmer = sequence.substr(sequence.size()-km, km);
            uint64_t hash_foward_kmer = hash_string(km, last_kmer, false);

            if (mode == "index"){
                if (right_seq.find(name) != right_seq.end()){
                    central_kmers_needed.insert(hash_foward_kmer);
                }
                else{
                    full_kmers_needed.insert(hash_foward_kmer);
                }
            }
            else
            {
                if (kmers.find(hash_foward_kmer) != kmers.end()){
                    //retrieve the full and central sequence from the kmers file
                    
                    string end_of_seq;
                    if (right_seq.find(name) != right_seq.end() && kmers[hash_foward_kmer].first != 1){ //if there was a right extension just take the central part
                        central_kmers_input.seekg(kmers[hash_foward_kmer].first);
                        std::getline(central_kmers_input, line2);
                        end_of_seq = line2;
                    }
                    else if (kmers[hash_foward_kmer].second != 1){
                        full_kmers_input.seekg(kmers[hash_foward_kmer].second);
                        std::getline(full_kmers_input, line2);
                        end_of_seq = line2;
                    }
                    //compute the overlap
                    int overlap = std::min(end_of_seq.size(), expanded_sequence.size());
                    int size_of_compared_seq = std::min(30,  overlap);
                    string exp_end = expanded_sequence.substr((int)expanded_sequence.size()-size_of_compared_seq, size_of_compared_seq);
                    while (overlap > 0 && end_of_seq.substr(overlap-size_of_compared_seq, size_of_compared_seq) != exp_end){
                        overlap--;
                        if (overlap < size_of_compared_seq){
                            exp_end = expanded_sequence.substr(expanded_sequence.size()-overlap, overlap);
                            size_of_compared_seq = overlap;
                        }
                    }
                    expanded_sequence += end_of_seq.substr(overlap, end_of_seq.size()-overlap);

                }
                else{
                    // cout << "WARNING (code 745) missing kmer " << last_kmer << "\n";
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
