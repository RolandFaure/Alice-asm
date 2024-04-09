// author : Roland Faure

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


#include "basic_graph_manipulation.h"
#include "robin_hood.h"
#include "clipp.h"

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

string version = "0.1.10";
string date = "2024-04-09";
string author = "Roland Faure";

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
void reduce(string input_file, string output_file, int context_length, int compression, int num_threads) {

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
                if (seq_num % 50000 == 0){
                    #pragma omp critical
                    {
                        //display the date and time
                        time_t now = time(0);
                        tm *ltm = localtime(&now);
                        cout << "[" << 1 + ltm->tm_mon << "/" << 1900 + ltm->tm_year << " " << 1 + ltm->tm_hour << ":" << 1 + ltm->tm_min << ":" << 1 + ltm->tm_sec << "]" << " Compressed " << seq_num << " reads" << endl;
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

                while (roll(hash_foward, hash_reverse, k, line, pos_end, pos_begin, true)){
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
void go_through_the_reads_again(string reads_file, string assemblyFile, int context_length, int compression, int km, unordered_map<string, pair<string,string>>& kmers, int num_threads){

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

        while (std::getline(input, line)){

            if (line[0] == '>')
            {
                if (seq_num % 50000 == 0){
                    #pragma omp critical
                    {
                        //display the date and time
                        time_t now = time(0);
                        tm *ltm = localtime(&now);
                        cout << "[" << 1 + ltm->tm_mon << "/" << 1900 + ltm->tm_year << " " << 1 + ltm->tm_hour << ":" << 1 + ltm->tm_min << ":" << 1 + ltm->tm_sec << "]" << " Processed " << seq_num << " reads" << endl;
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
                while (roll(hash_foward, hash_reverse, k, line, pos_end, pos_begin, pos_middle, true)){
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

void expand(string asm_reduced, string output, int km, int length_of_overlaps, unordered_map<string, pair<string,string>>& kmers){

    ifstream input(asm_reduced);

    ofstream out(output);
    out << "H\tVN:Z:1.0\n";

    int length_of_central_kmers = 1;
    int length_of_full_kmers = 1;

    //first index the first and last 10 bases of all contigs
    unordered_map<std::string, std::string> first_10;
    unordered_map<std::string, std::string> last_10;
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
            if (sequence.size() > 10 + length_of_overlaps){
                first_10[name] = sequence.substr(length_of_overlaps, std::min(10, (int)sequence.size()));
                last_10[name] = sequence.substr(std::max(0, (int)sequence.size()-10-length_of_overlaps), 10);
            }
            
        }
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

            if (orientation1 == "+"){                
                if (orientation2 == "+"){
                    right_seq[contig1] = first_10[contig2];
                    left_seq[contig2] = last_10[contig1];
                }
                else{
                    right_seq[contig1] = reverse_complement(last_10[contig2]);
                    right_seq[contig2] = reverse_complement(last_10[contig1]);
                }
            }
            else{
                if (orientation2 == "+"){
                    left_seq[contig1] = reverse_complement(first_10[contig2]);
                    left_seq[contig2] = reverse_complement(first_10[contig1]);
                }
                else{
                    left_seq[contig1] = last_10[contig2];
                    right_seq[contig2] = first_10[contig1];
                }
            }
        }
    }
    input.close();

    input.open(asm_reduced);
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

int main(int argc, char** argv)
{

    //use clipp to parse the command line
    bool help = false;
    string input_file, output_folder;
    string path_to_bcalm = "bcalm";
    string path_to_minimap = "minimap2";
    bool rescue = false;
    int min_abundance = 5;
    int order = 201;
    int compression = 20;
    int num_threads = 1;
    auto cli = (
        clipp::required("-r", "--reads").doc("input file (fasta/q)") & clipp::opt_value("r",input_file),
        clipp::required("-o", "--output").doc("output folder") & clipp::opt_value("o",output_folder),
        clipp::option("-t", "--threads").doc("number of threads [1]") & clipp::opt_value("t", num_threads),
        clipp::option("-m", "--min-abundance").doc("minimum abundance of kmer to consider solid [5]") & clipp::opt_value("m", min_abundance),
        clipp::option("-l", "--order").doc("order of MSR compression (odd) [201]") & clipp::opt_value("o", order),
        clipp::option("-c", "--compression").doc("compression factor [20]") & clipp::opt_value("c", compression),
        clipp::option("--bcalm").doc("path to bcalm [bcalm]") & clipp::opt_value("b", path_to_bcalm),
        clipp::option("--minimap2").doc("path to minimap2 [minimap2]") & clipp::opt_value("m", path_to_minimap),
        clipp::option("-v", "--version").call([]{ std::cout << "version " << version << "\nLast update: " << date << "\nAuthor: " << author << std::endl; exit(0); }).doc("print version and exit")

    );

    //ascii art of a cake:
//      __.--.__
//    .'  &  . '.\             
//   |    &     & |        
//   \______&____.'     
//    |.|.|.|.|.|                                        
//    `-|.|.|.|-' 

    //ascii art of a the Alice Assembler:
    //   _______ _                     _ _                                            _     _           
    //  |__   __| |              /\   | (_)              /\                          | |   | |          
    //     | |  | |__   ___     /  \  | |_  ___ ___     /  \   ___ ___  ___ _ __ ___ | |__ | | ___ _ __ 
    //     | |  | '_ \ / _ \   / /\ \ | | |/ __/ _ \   / /\ \ / __/ __|/ _ \ '_ ` _ \| '_ \| |/ _ \ '__|
    //     | |  | | | |  __/  / ____ \| | | (_|  __/  / ____ \\__ \__ \  __/ | | | | | |_) | |  __/ |   
    //     |_|  |_| |_|\___| /_/    \_\_|_|\___\___| /_/    \_\___/___/\___|_| |_| |_|_.__/|_|\___|_|   

    //ascii art of a bottle:
    //  (-)  
    //.-'-'-.
    //|-...-|
    //|;:.._|
    //| ... |
    //`-...-'

    //The Alice Assembler ascii art


    //Show the bottle left, then the cake right
    cout << "  (-)        _______ _                     _ _                                            _     _              " << "   __.--.__" << endl;
    cout << ".-'-'-.     |__   __| |              /\\   | (_)              /\\                          | |   | |             " << " .'  &  . '.\\" << endl;
    cout << "|-...-|        | |  | |__   ___     /  \\  | |_  ___ ___     /  \\   ___ ___  ___ _ __ ___ | |__ | | ___ _ __    " << "|    &    &  |" << endl;
    cout << "|;:.._|        | |  | '_ \\ / _ \\   / /\\ \\ | | |/ __/ _ \\   / /\\ \\ / __/ __|/ _ \\ '_ ` _ \\| '_ \\| |/ _ \\ '__|   "<< "\\______&____.'" << endl;
    cout << "| ... |        | |  | | | |  __/  / ____ \\| | | (_|  __/  / ____ \\\\__ \\__ \\  __/ | | | | | |_) | |  __/ |      " << " |.|.|.|.|.| " << endl;
    cout << "`-...-'        |_|  |_| |_|\\___| /_/    \\_\\_|_|\\___\\___| /_/    \\_\\___/___/\\___|_| |_| |_|_.__/|_|\\___|_|      " << " `-|.|.|.|-'"<< endl;
    cout << endl;
     
    cout << "Command line: ";
    for (int i = 0 ; i < argc ; i++){
        cout << argv[i] << " ";
    }
    cout << endl;

    cout << "Alice Assembler version " << version << "\nLast update: " << date << "\nAuthor: " << author << endl << endl;


    if(!clipp::parse(argc, argv, cli)) {
        cout << "Could not parse the arguments" << endl;
        cout << clipp::make_man_page(cli, argv[0]);
        exit(1);
    }

    if (order % 2 == 0){
        cerr << "ERROR: order must be odd\n";
        exit(1);
    }
    int context_length = (order-1)/2;

    //make sure the output folder ends with a /
    if (output_folder[output_folder.size()-1] != '/'){
        output_folder += "/";
    }
    string tmp_folder = output_folder + "tmp/";

    //record time now to measure the time of the whole process
    auto start = std::chrono::high_resolution_clock::now();

    //check if the output folder exists
    string command = "mkdir -p " + output_folder + " 2> /dev/null";
    system(command.c_str());

    //create the tmp folder
    command = "mkdir -p " + tmp_folder + " 2> /dev/null";
    system(command.c_str());
    

    string path_src = argv[0];
    path_src = path_src.substr(0, path_src.find_last_of("/")); //strip the /reduce
    path_src = path_src.substr(0, path_src.find_last_of("/")); //strip the /build

    //check if the path to bcalm is correct
    string command_bcalm = path_to_bcalm + " --help 2> /dev/null > /dev/null";
    auto bcalm_ok = system(command_bcalm.c_str());
    if (bcalm_ok != 0){
        cerr << "ERROR: bcalm not found using command line " << command_bcalm << "\n";
        cerr << "Please specify a valid path using --bcalm\n" << endl;
        exit(1);
    }
    //check if the path to minimap is correct
    command = path_to_minimap + " --help 2> /dev/null > /dev/null";
    auto ok = system(command.c_str());
    if (ok != 0){
        cerr << "ERROR: minimap2 not found using command line " << command << "\n";
        cerr << "Please specify a valid path using --minimap2\n" << endl;
        exit(1);
    }

    std::string path_convertToGFA = "python " + path_src + "/bcalm/scripts/convertToGFA.py";

    int km = 31; //size of the kmer used to do the expansion. Must be >21
    vector<int> values_of_k = {16,31,41,71}; //size of the kmer used to build the graph (min >= km)

    //if the input file is a fastq file, convert it to fasta
    if (input_file.substr(input_file.find_last_of('.')+1) == "fastq" || input_file.substr(input_file.find_last_of('.')+1) == "fq"){
        string fasta_file = input_file.substr(0, input_file.find_last_of('.')) + ".fasta";
        string command = "sed -n '1~4s/^@/>/p;2~4p' " + input_file + " > " + tmp_folder + fasta_file;
        auto ok = system(command.c_str());
        if (ok != 0){
            cerr << "ERROR: fastq_to_fasta failed\n";
            exit(1);
        }
        input_file = tmp_folder+fasta_file;
    }

    string compressed_file = tmp_folder+"compressed.fa";
    string merged_gfa = tmp_folder+"bcalm.unitigs.shaved.merged.gfa";

    cout << "==== Step 1: MSR compression of the reads ====\n";
    
    reduce(input_file, compressed_file, context_length, compression, num_threads);

    cout << "Done compressing reads, the compressed reads are in " << compressed_file << "\n" << endl;

    bool hifiasm = false;
    string untangled_gfa;
    if (!hifiasm){

        cout << "==== Step 2: Iterative DBG assemby of the compressed reads with increasing k ====\n";
        //iterative k to have good contiguity
        for (auto kmer_len: values_of_k){
            // launch bcalm        
            cout << " - Launching assembly with k=" << kmer_len << endl;
            cout << "    - Unitig generation with bcalm" << endl;
            string bcalm_command = path_to_bcalm + " -in " + compressed_file + " -kmer-size "+std::to_string(kmer_len)+" -abundance-min " 
                + std::to_string(min_abundance) + " -out "+tmp_folder+"bcalm > "+tmp_folder+"bcalm.log 2>&1";
            auto bcalm_ok = system(bcalm_command.c_str());
            if (bcalm_ok != 0){
                cerr << "ERROR: bcalm failed\n";
                cout << bcalm_command << endl;
                exit(1);
            }

            // convert to gfa
            cout << "    - Converting result to GFA" << endl;
            string unitig_file_fa = tmp_folder+"bcalm.unitigs.fa";
            string unitig_file_gfa = tmp_folder+"bcalm.unitigs.gfa";
            string convert_command = path_convertToGFA + " " + unitig_file_fa + " " + unitig_file_gfa +" "+ std::to_string(kmer_len) + " > " + tmp_folder + "convertToGFA.log 2>&1";
            system(convert_command.c_str());

            // shave the resulting graph
            cout << "    - Shaving the graph of small dead ends" << endl;
            string shaved_gfa = tmp_folder+"bcalm.unitigs.shaved.gfa";
            shave(unitig_file_gfa, shaved_gfa, 2*kmer_len-1);

            //merge the adjacent contigs
            cout << "    - Merging resulting contigs" << endl;
            merge_adjacent_contigs_BCALM(shaved_gfa, merged_gfa, kmer_len, path_to_bcalm, path_convertToGFA, tmp_folder);

            //take the contigs of bcalm.unitigs.shaved.merged.unzipped.gfa and put them in a fasta file min_abundance times, and concatenate with compressed_file
            cout << "    - Concatenating the contigs to the reads to relaunch assembly with higher k" << endl;
            
            //open both compressed_file and bcalm.unitigs.shaved.merged.unzipped.gfa
            ofstream input_compressed(compressed_file, std::ios_base::app);
            ifstream input_graph(merged_gfa);
            string line;
            while (std::getline(input_graph, line))
            {
                if (line[0] == 'S')
                {
                    string name;
                    string dont_care;
                    string sequence;
                    std::stringstream ss(line);
                    ss >> dont_care >> name >> sequence;
                    for (int i = 0 ; i < min_abundance ; i++){
                        input_compressed << ">" << name << "\n";
                        input_compressed << sequence << "\n";
                    }
                }
            }
            input_compressed.close();
            input_graph.close();
        }

        cout << "Done with the iterative assembly, the final graph is in " << merged_gfa << "\n" << endl;

        cout << "==== Step 3: Untangling the final compressed assembly ====\n";

        //sort the gfa to have S lines before L lines
        cout << " - Sorting the GFA" << endl;
        sort_GFA(merged_gfa);

        //untangle the graph to improve contiguity
        cout << " - Aligning the reads to the graph" << endl;
        string gaf_file = tmp_folder+"bcalm.unitigs.shaved.merged.unzipped.gaf";
        unordered_map<string,float> coverages;
        create_gaf_from_unitig_graph(merged_gfa, values_of_k[values_of_k.size()-1], compressed_file, gaf_file, coverages);
        add_coverages_to_graph(merged_gfa, coverages);
        
        cout << " - Untangling the graph with GraphUnzip" << endl;
        untangled_gfa = tmp_folder+"bcalm.unitigs.shaved.merged.unzipped.untangled.gfa";
        string command_unzip = "python " + path_src + "/GraphUnzip/graphunzip.py unzip -R -l " + gaf_file + " -g " + merged_gfa + " -o " + untangled_gfa + " > " + tmp_folder + "graphunzip.log 2>&1";
        auto unzip_ok = system(command_unzip.c_str());
        if (unzip_ok != 0){
            cerr << "ERROR: unzip failed\n";
            exit(1);
        }
        cout << "Done untangling the graph, the final compressed graph is in " << untangled_gfa << "\n" << endl;
    }
    else{
        //assemble the compressed reads with hifiasm
        cout << "==== Step 2: Assembly of the compressed reads with hifiasm ====\n";
        string hifiasm_output = tmp_folder;
        string command_hifiasm = "hifiasm -o " + hifiasm_output + "hifiasm -t " + std::to_string(num_threads) + " " + compressed_file + " > " + tmp_folder + "hifiasm.log 2>&1";
        untangled_gfa = hifiasm_output + "hifiasm.bp.p_ctg.gfa";
        auto hifiasm_ok = system(command_hifiasm.c_str());
        if (hifiasm_ok != 0){
            cerr << "ERROR: hifiasm failed after running command line\n";
            cerr << command_hifiasm << endl;
            exit(1);
        }
    }

    cout << "==== Step 4: Inflating back the assembly to non-compressed space ====\n";

    //now let's parse the gfa file and decompress it
    cout << " - Parsing the reads to map compressed kmers with uncompressed sequences\n";
    //time the next function
    unordered_map<string, pair<string,string>> kmers;
    go_through_the_reads_again(input_file, untangled_gfa, context_length, compression, km, kmers, num_threads);
    string decompressed_file = tmp_folder+"bcalm.unitigs.shaved.merged.unzipped.decompressed.gfa";
    cout << " - Reconstructing the uncompressed assembly" << endl;
    expand(untangled_gfa, decompressed_file, km, values_of_k[values_of_k.size()-1]-1, kmers);

    //test pop_and_shave_homopolymer_errors //not necessary now that the compression is done with HPC
    // string pop_and_shave_file = "bcalm.unitigs.shaved.merged.unzipped.decompressed.pop_and_shaved.gfa";
    // remove_homopolymer_errors(decompressed_file, pop_and_shave_file, path_to_minimap);
    string pop_and_shave_file = decompressed_file;
    // exit(0);

    //test compute_eact_cigar
    string output_file = output_folder + "assembly.gfa";
    cout << " - Computing the exact overlaps between the contigs\n";
    compute_exact_CIGARs(pop_and_shave_file, output_file, 2000, (values_of_k[values_of_k.size()-1]-1)*compression);

    //convert to fasta
    
    cout << " - Converting the assembly to fasta\n";
    gfa_to_fasta(output_file, output_file.substr(0, output_file.find_last_of('.')) + ".fasta");

    //remove the temporary file
    // command = "cd "+ tmp_folder +" && ls | grep -v 'shaved' | xargs rm && cd ../../";
    // system(command.c_str());


    cout << "\nDone, the final assembly is in " << output_file << "\n" << endl;
    cout << "Total time: " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count() << "s\n";

    return 0;
}

