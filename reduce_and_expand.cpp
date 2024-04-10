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
 * @param context_length 
 * @param compression 
 * @param km 
 * @param min_abundance 
 * @param kmers maps a kmer to the uncompressed seq
 **/
void reduce(string input_file, string output_file, int context_length, int compression, int num_threads, bool homopolymer_compression) {

    int k = 2*context_length + 1;

    std::ifstream input(input_file,std::ios::binary | std::ios::ate);
    if (!input.is_open())
    {
        std::cout << "Could not open file " << input_file << std::endl;
        exit(1);
    }

    std::streamoff file_size = input.tellg();
    unsigned long long size_of_chunk = 100000000;
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

        while (std::getline(input, line))
        {
            if (line[0] == '>')
            {
                out << line << "\n";
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
            }
            else if (next_line_is_seq){
                //let's launch the foward and reverse rolling hash
                uint64_t hash_foward = 0;
                uint64_t hash_reverse = 0;
                size_t pos_end = 0;
                long pos_begin = -k;

                while (roll(hash_foward, hash_reverse, k, line, pos_end, pos_begin, homopolymer_compression)){
                    if (pos_begin>=0){
                        if (hash_foward<hash_reverse && hash_foward % compression == 0){
                            out << "ACGT"[(hash_foward/compression)%4];
                        }
                        else if (hash_foward>=hash_reverse && hash_reverse % compression == 0){
                            out << "TGCA"[(hash_reverse/compression)%4];
                        }
                    }
                }
                out << "\n";
            }
            //if we went beyond chunk+1 * size_of_chunk, we can stop
            if (input.tellg() > (chunk+1)*size_of_chunk){
                break;
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
 * @param km size of the kmer used to do the expansion
 * @param kmers Associates a compressed sequence to its central uncompressed part and its full uncompressed part
 */
void go_through_the_reads_again(string reads_file, string assemblyFile, int context_length, int compression, int km, unordered_map<string, pair<string,string>>& kmers, int num_threads, bool homopolymer_compression){

    unordered_set<string> kmers_in_assembly;

    ifstream input_asm(assemblyFile);
    string line;
    while (std::getline(input_asm, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            for (int i = 0 ; i <= sequence.size()-km ; i++){
                string kmer = sequence.substr(i, km);
                kmers_in_assembly.insert(kmer);
            }
        }
    }
    input_asm.close();

    unordered_map<string, unordered_map<string, int>> kmer_count;
    unordered_map<string, bool> confirmed_kmers;

    int k = 2*context_length + 1;

    ifstream input2(reads_file, std::ios::binary | std::ios::ate);
    std::streamoff file_size = input2.tellg();
    unsigned long long int size_of_chunk = 100000000;
    input2.close();

    //go through the fasta file and compute the hash of all the kmer using homecoded ntHash
    int seq_num = 0;
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
    for (int chunk = 0 ; chunk <= file_size/size_of_chunk ; chunk++){

        std::ifstream input(reads_file);
        input.seekg(chunk*size_of_chunk);
    
        std::string line;
        bool next_line_is_seq = false;
        string read_name;
        long long output_limit = 0;

        while (std::getline(input, line)){

            if (line[0] == '>')
            {
                if (seq_num > output_limit && omp_get_thread_num() == 0){
                    #pragma omp critical
                    {
                        //display the date and time
                        time_t now = time(0);
                        tm *ltm = localtime(&now);
                        cout << "[" << 1 + ltm->tm_mon << "/" << 1900 + ltm->tm_year << " " << ltm->tm_hour << ":" << ltm->tm_min << ":" << ltm->tm_sec << "]" << " Processed " << seq_num << " reads" << endl;
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
                long pos_begin = -k;
                long pos_middle = -(k+1)/2;

                vector<int> positions_sampled (0);
                std::string kmer = string(km, 'N');
                string rkmer;
                while (roll(hash_foward, hash_reverse, k, line, pos_end, pos_begin, pos_middle, homopolymer_compression)){
                    if (pos_begin>=0){

                        if ((hash_foward<hash_reverse && hash_foward % compression == 0) || (hash_foward>=hash_reverse && hash_reverse % compression == 0)){

                            if (hash_foward<hash_reverse){
                                kmer = kmer.substr(1,kmer.size()-1) + "ACGT"[(hash_foward/compression)%4];
                            }
                            else{
                                kmer = kmer.substr(1,kmer.size()-1) + "TGCA"[(hash_reverse/compression)%4];
                            }

                            rkmer = reverse_complement(kmer); 
                            positions_sampled.push_back(pos_middle);

                            // if (kmer == "GTACCACTGACCTTACATATCGATTGTTTAA"){
                            //     cout << "FOUND " << line.substr(pos_middle -1, 5) << " in " << read_name << endl;
                            //     cout << positions_sampled[positions_sampled.size()-km+10] << " " << positions_sampled[positions_sampled.size()-1-10] << endl;
                            //     cout << line.substr(positions_sampled[positions_sampled.size()-km+10], positions_sampled[positions_sampled.size()-1-10] - positions_sampled[positions_sampled.size()-km+10]+1) << endl;
                            //     cout << "confirmed " << confirmed_kmers["GTACCACTGACCTTACATATCGATTGTTTAA"] << " " << kmer_count["GTACCACTGACCTTACATATCGATTGTTTAA"] << endl;
                            //     if (confirmed_kmers["GTACCACTGACCTTACATATCGATTGTTTAA"]){
                            //         cout << "confirmed " << kmers["GTACCACTGACCTTACATATCGATTGTTTAA"].first << endl;
                            //         exit(1);
                            //     }

                            // }
                            // else if (rkmer == "GCCATGACAACCTCTCGCCTTCTGGAGCCGT"){
                            //     cout << "FOUNDdd " << line << endl;
                            //     exit(1);
                            // }

                            if (positions_sampled.size() >= km && (kmers_in_assembly.find(kmer) != kmers_in_assembly.end() || kmers_in_assembly.find(rkmer) != kmers_in_assembly.end())){
                                
                                string canonical_kmer = min(kmer, rkmer);

                                #pragma omp critical
                                {
                                    if (confirmed_kmers[canonical_kmer] == false){

                                        if (kmer_count.find(canonical_kmer) == kmer_count.end()){
                                            kmer_count[canonical_kmer] = {};
                                            kmers[kmer] = {"",""}; //first member is the central, "sure" part, the second is the full sequence, but potentially with a little noise at the ends
                                            kmers[rkmer] = {"",""};
                                            confirmed_kmers[canonical_kmer] = false;
                                        }
                                        
                                        string central_seq = line.substr(positions_sampled[positions_sampled.size()-km+10], positions_sampled[positions_sampled.size()-1-10] - positions_sampled[positions_sampled.size()-km+10]+1);
                                        string reverse_central_seq = reverse_complement(central_seq);
                                        string canonical_central_seq = min(central_seq, reverse_central_seq);
                                        if (kmer_count[canonical_kmer].find(canonical_central_seq) == kmer_count[canonical_kmer].end()){
                                            kmer_count[canonical_kmer][canonical_central_seq] = 0;
                                        }
                                        kmer_count[canonical_kmer][canonical_central_seq]++;
                                        
                                        if (kmer_count[canonical_kmer][canonical_central_seq] > 3){
                                            string full_seq = line.substr(positions_sampled[positions_sampled.size()-km], positions_sampled[positions_sampled.size()-1] - positions_sampled[positions_sampled.size()-km]+1);
                                            kmers[kmer] = {central_seq, full_seq};
                                            kmers[rkmer] = {reverse_complement(central_seq), reverse_complement(full_seq)};
                                            confirmed_kmers[canonical_kmer] = true;
                                            kmer_count[canonical_kmer].clear();
                                        }
                                    }
                                }
                            }
                        }
                    }
                    
                }
                //if we went beyond chunk+1 * size_of_chunk, we can stop
                if (input.tellg() > (chunk+1)*size_of_chunk){
                    break;
                }
            }
        }
        input.close();

    }
}

/**
 * @brief Takes in the reduced assembly and expands it
 * 
 * @param asm_reduced 
 * @param output 
 * @param km size of k used for compression
 * @param length_of_overlaps size of k used for assembly -1 
 * @param kmers dictionnary mapping compressed kmers to their uncompressed central part and full part
 */
void expand(string asm_reduced, string output, int km, int length_of_overlaps, unordered_map<string, pair<string,string>>& kmers){

    ifstream input(asm_reduced);

    ofstream out(output);
    out << "H\tVN:Z:1.0\n";

    int length_of_central_kmers = 1;
    int length_of_full_kmers = 1;

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
            // if (sequence.size() > 10 + length_of_overlaps){
            //     first_10[name] = sequence.substr(length_of_overlaps, std::min(10, (int)sequence.size()));
            //     last_10[name] = sequence.substr(std::max(0, (int)sequence.size()-10-length_of_overlaps), 10);
            // }
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

    input.open(asm_reduced);
    long long position_in_file;
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
            
            //expand the sequence (focusing on the central part of each kmer and thus missing the two ends)
            string expanded_sequence = "";
            int i = 0;
            for (i = 0; i <= sequence.size()-km; i+= km-20-1){ //-20 because we only take the central part of each kmer
                string kmer = sequence.substr(i, km);
                if (kmers.find(kmer) != kmers.end()){
                    length_of_central_kmers = kmers[kmer].first.size();
                    length_of_full_kmers = kmers[kmer].second.size();
                    if (i == 0 || kmers[kmer].first.size() == 0){
                        expanded_sequence += kmers[kmer].first;
                    }
                    else{ //don't append the first base, it has already been appended in the previous kmer
                        expanded_sequence.append(kmers[kmer].first.begin()+1, kmers[kmer].first.end());
                    }
                }
                else{
                    cerr << "WARNING (code 554) not found kmer " << kmer << "\n";
                    // cerr << line << endl;
                    // exit(1);
                    //create a sequence of Ns of sze length_of_central_kmers
                    expanded_sequence += string(length_of_central_kmers, 'N');
                }
            }
            
            // create the beginning of the sequence if there was no left extension
            if (left_seq.find(name) == left_seq.end()){
                string first_kmer = sequence.substr(0, km);
                if (kmers.find(first_kmer) != kmers.end()){
                    string beginning_of_seq = kmers[first_kmer].second;

                    //compute the overlap
                    int overlap = beginning_of_seq.size();
                    string exp_start = expanded_sequence.substr(0, std::min(30, overlap));
                    while (overlap > 0 && beginning_of_seq.substr(beginning_of_seq.size()-overlap, exp_start.size()) != exp_start){
                        overlap--;
                        if (overlap < 30){
                            exp_start = expanded_sequence.substr(0, overlap);
                        }
                    }
                    expanded_sequence = beginning_of_seq.substr(0, beginning_of_seq.size()-overlap) + expanded_sequence;
                }
                else{
                    cerr << "WARNING (code 555) not found kmer " << first_kmer << "\n";
                }
            }

            // finish the sequence
            string last_kmer = sequence.substr(sequence.size()-km, km);
            if (kmers.find(last_kmer) != kmers.end()){
                string end_of_seq = kmers[last_kmer].second;
                if (right_seq.find(name) != right_seq.end()){ //if there was a right extension just take the central part
                    end_of_seq = kmers[last_kmer].first;
                }
                //compute the overlap
                int overlap = end_of_seq.size();

                string exp_end = expanded_sequence.substr(std::max((int)expanded_sequence.size()-30,0), std::min((int)expanded_sequence.size(), 30));
                while (overlap > 0 && end_of_seq.substr(overlap-exp_end.size(), exp_end.size()) != exp_end){
                    overlap--;
                    if (overlap < 30){
                        exp_end = expanded_sequence.substr(expanded_sequence.size()-overlap, overlap);
                    }
                }

                expanded_sequence += end_of_seq.substr(overlap, end_of_seq.size()-overlap);
            }
            else{
                cerr << "WARNING (code 557) not found kmer " << last_kmer << "\n";
            }
 
            out << "S\t" << name << "\t" << expanded_sequence;
            while (ss >> sequence){
                out << "\t" << sequence;
            }
            out << "\n";
        }
        else{
            out << line << "\n";
        }
    }
}

