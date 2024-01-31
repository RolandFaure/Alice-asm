// author : Roland Faure

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <unordered_set>

#include <nthash/nthash.hpp>

#include "basic_graph_manipulation.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::unordered_map;
using std::pair;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::unordered_set;

string version = "0.1.1";
string date = "2024-01-31";
string author = "Roland Faure";

/**
 * @brief 
 * 
 * @param input_file 
 * @param output_file 
 * @param context_length 
 * @param compression 
 * @param km 
 * @param min_abundance 
 * @param kmers maps a kmer to the uncompressed seq
 **/
void reduce(string input_file, string output_file, int context_length, int compression, int km, unordered_map<string, pair<string,string> > &kmers) {

    std::ifstream input(input_file);
    if (!input.is_open())
    {
        std::cout << "Could not open file " << input_file << std::endl;
        exit(1);
    }
    unordered_map<string, int> kmer_count;
    unordered_map<string, bool> confirmed_kmers;

    std::ofstream out(output_file);

    int k = 2*context_length + 1;
    string seed = string(context_length, '1') + "0" + string(context_length, '1');
    vector<string> seeds = {seed};

    //go through the fasta file and compute the hash of all the kmer using ntHash
    std::string line;

    long identical = 0;
    long different = 0;

    int seq_num = 0;
    while (std::getline(input, line))
    {
        if (line[0] == '>')
        {
            out << line << "\n";
            if (seq_num % 100 == 0){
                cout << "Compressing read " << seq_num << "\r";
            }
            seq_num++;
        }
        else{

            nthash::SeedNtHash nth(line, seeds, 1, k);
            vector<int> positions_sampled (0);
            int pos = 0;
            std::string kmer = string(km, 'N');
            string rkmer;
            while (nth.roll()) {        
                // cout << "hash of " << line.substr(pos, k) << " is " << nth.hashes()[0] << "\n";
                if ( pos < line.size() - k && nth.hashes()[0] % compression == 0){
                    out << line[pos+context_length];
                }
                pos++;
                if (pos+k > line.size()){
                    break; //or else it will roll to the beginning of the sequence
                }
            }
            out << "\n";
        }
    }
    input.close();
    out.close();
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
void go_through_the_reads_again(string reads_file, string assemblyFile, int context_length, int compression, int km, unordered_map<string, pair<string,string>>& kmers){

    unordered_set<string> kmers_in_assembly;

    ifstream input(assemblyFile);
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
            for (int i = 0 ; i <= sequence.size()-km ; i++){
                string kmer = sequence.substr(i, km);
                kmers_in_assembly.insert(kmer);
            }
        }
    }
    input.close();

    input.open(reads_file);
    unordered_map<string, int> kmer_count;
    unordered_map<string, bool> confirmed_kmers;

    int k = 2*context_length + 1;
    string seed = string(context_length, '1') + "0" + string(context_length, '1');
    vector<string> seeds = {seed};

    //go through the fasta file and compute the hash of all the kmer using ntHash
    int seq_num = 0;
    while (std::getline(input, line))
    {
        if (line[0] == '>')
        {
            if (seq_num % 100 == 0){
                cout << "Compressing read " << seq_num << "\r";
                // cout << "nanaaaammmmma " << line << endl;
            }
            seq_num++;
        }
        else{

            nthash::SeedNtHash nth(line, seeds, 1, k);
            vector<int> positions_sampled (0);
            int pos = 0;
            std::string kmer = string(km, 'N');
            string rkmer;
            while (nth.roll()) {        
                if ( pos < line.size() - k && nth.hashes()[0] % compression == 0){
                    // cout << "hash of " << line.substr(pos, k) << " is " << nth.hashes()[0] << "\n";
                    // exit(1);
                    positions_sampled.push_back(pos+context_length);
                    kmer = kmer.substr(1,kmer.size()-1) + line[pos+context_length];
                    rkmer = reverse_complement(kmer); 

                    if (positions_sampled.size() >= km && (kmers_in_assembly.find(kmer) != kmers_in_assembly.end() || kmers_in_assembly.find(rkmer) != kmers_in_assembly.end())){
                        
                        string canonical_kmer = min(kmer, rkmer);

                        if (kmer_count.find(canonical_kmer) == kmer_count.end()){
                            kmer_count[canonical_kmer] = 0;
                            kmers[kmer] = {"",""}; //first member is the central, "sure" part, the second is the full sequence, but potentially with a little noise at the ends
                            kmers[rkmer] = {"",""};
                        }

                        kmer_count[canonical_kmer]++;
                        if (kmer_count[canonical_kmer] == 1){
                            string central_seq = line.substr(positions_sampled[positions_sampled.size()-km+10], positions_sampled[positions_sampled.size()-1-10] - positions_sampled[positions_sampled.size()-km+10]);
                            string full_seq = line.substr(positions_sampled[positions_sampled.size()-km], positions_sampled[positions_sampled.size()-1] - positions_sampled[positions_sampled.size()-km]);
                            kmers[kmer] = {central_seq, full_seq};
                            //since we have to exclude the last base, the reverse complement is slightly different from the foward
                            central_seq = line.substr(positions_sampled[positions_sampled.size()-km+10]+1, positions_sampled[positions_sampled.size()-1-10] - positions_sampled[positions_sampled.size()-km+10]);
                            full_seq = line.substr(positions_sampled[positions_sampled.size()-km]+1, positions_sampled[positions_sampled.size()-1] - positions_sampled[positions_sampled.size()-km]);
                            kmers[rkmer] = {reverse_complement(central_seq), reverse_complement(full_seq)};
                        }
                        else if (kmer_count[canonical_kmer] > 1){

                            if (!confirmed_kmers[canonical_kmer]){
                                string central_seq = line.substr(positions_sampled[positions_sampled.size()-km+10], positions_sampled[positions_sampled.size()-1-10] - positions_sampled[positions_sampled.size()-km+10]);
                                if (kmers[kmer].first != central_seq){
                                    string full_seq = line.substr(positions_sampled[positions_sampled.size()-km], positions_sampled[positions_sampled.size()-1] - positions_sampled[positions_sampled.size()-km]);
                                    kmers[kmer].first = central_seq;
                                    kmers[kmer].second = full_seq;
                                    central_seq = line.substr(positions_sampled[positions_sampled.size()-km+10]+1, positions_sampled[positions_sampled.size()-1-10] - positions_sampled[positions_sampled.size()-km+10]);
                                    full_seq = line.substr(positions_sampled[positions_sampled.size()-km]+1, positions_sampled[positions_sampled.size()-1] - positions_sampled[positions_sampled.size()-km]);
                                    kmers[rkmer].first = reverse_complement(central_seq);
                                    kmers[rkmer].second = reverse_complement(full_seq);
                                }
                                else{
                                    confirmed_kmers[canonical_kmer] = true;
                                }
                            }
                        }
                    }
                }                
                
                pos++;
                if (pos+k > line.size()){
                    break; //or else it will roll to the beginning of the sequence
                }
            }
        }
    }
}

void expand(string asm_reduced, string output, int km, unordered_map<string, pair<string,string>>& kmers){

    ifstream input(asm_reduced);

    ofstream out(output);
    out << "H\tVN:Z:1.0\n";

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
            
            //expand the sequence (focusing on the central part of each kmer and thus missing the two ends)
            string expanded_sequence = "";
            int i = 0;
            for (i = 0; i <= sequence.size()-km; i+= km-20-1){ //-20 because we only take the central part of each kmer
                string kmer = sequence.substr(i, km);
                if (name == "240--1"){
                    cout << "expanding 262 with " << kmer << " " << kmers[kmer].first << "\n";
                }
                if (name == "114--1"){
                    cout << "expanding 100 with " << kmer << " " << kmers[kmer].first << "\n";
                }
                if (kmers.find(kmer) != kmers.end()){
                    expanded_sequence += kmers[kmer].first;
                }
                else{
                    cerr << "ERROR 1 not found kmer " << kmer << "\n";
                    cerr << line << endl;
                    exit(1);
                }
            }
            //begin the sequence
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
                if (name == "240--1"){
                    cout << "beginning 262 with " << first_kmer << " " << beginning_of_seq << " " << beginning_of_seq.substr(0, beginning_of_seq.size()-overlap) << "\n";
                }
                if (name == "114--1"){
                    cout << "beginning 100 with " << first_kmer << " " << beginning_of_seq << " " << beginning_of_seq.substr(0, beginning_of_seq.size()-overlap) << "\n";
                }
                expanded_sequence = beginning_of_seq.substr(0, beginning_of_seq.size()-overlap) + expanded_sequence;
            }
            else{
                cerr << "ERROR 2 not found kmer " << first_kmer << "\n";
                cerr << line << endl;
                exit(1);
            }
            // finish the sequence
            string last_kmer = sequence.substr(sequence.size()-km, km);
            if (kmers.find(last_kmer) != kmers.end()){
                string end_of_seq = kmers[last_kmer].second;
                //compute the overlap
                int overlap = end_of_seq.size();
                string exp_end = expanded_sequence.substr(std::min((int)expanded_sequence.size()-30,0), std::min((int)expanded_sequence.size(), 30));
                while (overlap > 0 && end_of_seq.substr(overlap-exp_end.size(), exp_end.size()) != exp_end){
                    overlap--;
                    if (overlap < 30){
                        exp_end = expanded_sequence.substr(expanded_sequence.size()-overlap, overlap);
                    }
                }
                if (name == "240--1"){
                    cout << "finishing 262 with " << last_kmer << " " << end_of_seq << " " << end_of_seq.substr(overlap, end_of_seq.size()-overlap) << "\n";
                }
                if (name == "114--1"){
                    cout << "finishing 100 with " << last_kmer << " " << end_of_seq << " "<< end_of_seq.substr(overlap, end_of_seq.size()-overlap) << "\n";
                }
                expanded_sequence += end_of_seq.substr(overlap, end_of_seq.size()-overlap);
            }
            else{
                cerr << "ERROR 3 not found kmer " << last_kmer << "\n";
                cerr << line << endl;
                exit(1);
            }

            out << "S\t" << name << "\t" << expanded_sequence;
            string otherTags;
            while (ss >> otherTags){
                out << "\t" << otherTags;
            }
            out << "\n";
        }
        else{
            out << line << "\n";
        }
    }
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

int main(int argc, char** argv)
{

    if (argc != 5)
    {
        std::cout << "Usage: " << argv[0] << " <input_file> <output_file> <context_length> <compression>" << std::endl;
        return 1;
    }

    string path_src = argv[0];
    path_src = path_src.substr(0, path_src.find_last_of("/")); //strip the /reduce
    path_src = path_src.substr(0, path_src.find_last_of("/")); //strip the /build
    string input_file = argv[1];
    string output_file = argv[2];
    int context_length = atoi(argv[3]);
    int compression = atoi(argv[4]);
    int min_abundance = 10;
    int km = 71; //size of the kmer used to do the expansion
    vector<int> values_of_k = {17, 41, 71}; //size of the kmer used to build the graph (max >= km)

    //if the input file is a fastq file, convert it to fasta
    if (input_file.substr(input_file.find_last_of('.')+1) == "fastq" || input_file.substr(input_file.find_last_of('.')+1) == "fq"){
        string fasta_file = input_file.substr(0, input_file.find_last_of('.')) + ".fasta";
        string command = "sed -n '1~4s/^@/>/p;2~4p' " + input_file + " > " + fasta_file;
        auto ok = system(command.c_str());
        if (ok != 0){
            cerr << "ERROR: fastq_to_fasta failed\n";
            exit(1);
        }
        input_file = fasta_file;
    }

    string compressed_file = input_file + ".compressed";

    unordered_map<string, pair<string,string>> kmers;
    reduce(input_file, compressed_file, context_length, compression, km, kmers);
    cout << "finished reducing\n";

    //iterative k to have good contiguity
    for (auto kmer_len: values_of_k){
        // launch bcalm
        cout << "Launching bcalm with k=" << kmer_len << "\n";
        string bcalm_command = path_src + "/bcalm/build/bcalm -in " + compressed_file + " -kmer-size "+std::to_string(kmer_len)+" -abundance-min " + std::to_string(min_abundance) + " -out bcalm";
        auto bcalm_ok = system(bcalm_command.c_str());
        if (bcalm_ok != 0){
            cerr << "ERROR: bcalm failed\n";
            cout << bcalm_command << endl;
            exit(1);
        }

        // convert to gfa
        cout << "Launching convertToGFA\n";
        string convert_command = path_src + "/bcalm/scripts/convertToGFA.py bcalm.unitigs.fa bcalm.unitigs.gfa " + std::to_string(kmer_len);
        system(convert_command.c_str());

        // shave the resulting graph
        cout << "Launching shave\n";
        shave("bcalm.unitigs.gfa", "bcalm.unitigs.shaved.gfa", 2*km-1);

        //convert bcaml.unitigs.shaved.gfa to fasta
        cout << "Convert to fasta bcalm.unitigs.shaved.gfa\n";
        gfa_to_fasta("bcalm.unitigs.shaved.gfa", "bcalm.unitigs.shaved.fasta");

        //to merge, simply make a unitig graph from bcalm.unitigs.shaved.gfa and then convert it to gfa
        cout << "Creating shaved unitig graph\n";
        string command_unitig_graph =path_src + "/bcalm/build/bcalm -in bcalm.unitigs.shaved.fasta -kmer-size "+std::to_string(kmer_len)+" -abundance-min 1 -out bcalm.shaved.merged";
        auto unitig_graph_ok = system(command_unitig_graph.c_str());
        if (unitig_graph_ok != 0){
            cerr << "ERROR: unitig graph failed\n";
            cout << command_unitig_graph << endl;
            exit(1);
        }

        //convert to gfa
        cout << "Launching convertToGFA\n";
        string convert_command2 = path_src + "/bcalm/scripts/convertToGFA.py bcalm.shaved.merged.unitigs.fa bcalm.unitigs.shaved.merged.gfa " + std::to_string(kmer_len);
        system(convert_command2.c_str());

        //sort the gfa to have S lines before L lines
        cout << "Sorting GFA\n";
        sort_GFA("bcalm.unitigs.shaved.merged.gfa");

        //untangle the graph to improve contiguity
        cout << "Creating GAF file\n";
        string gaf_file = "bcalm.unitigs.shaved.merged.unzipped.gaf";
        create_gaf_from_unitig_graph("bcalm.unitigs.shaved.merged.gfa", kmer_len, compressed_file, gaf_file);
        
        cout << "Untangling GFA\n";
        string command_unzip = "python " + path_src + "/GraphUnzip/graphunzip.py unzip -R -l bcalm.unitigs.shaved.merged.unzipped.gaf -g bcalm.unitigs.shaved.merged.gfa -o bcalm.unitigs.shaved.merged.unzipped.gfa";
        cout << "command is " << command_unzip << endl;
        auto unzip_ok = system(command_unzip.c_str());
        if (unzip_ok != 0){
            cerr << "ERROR: unzip failed\n";
            exit(1);
        }

        //take the contigs of bcalm.unitigs.shaved.merged.unzipped.gfa and put them in a fasta file min_abundance times, and concatenate with compressed_file
        cout << "Appending to the compressed reads the contigs found\n";
        
        //open both compressed_file and bcalm.unitigs.shaved.merged.unzipped.gfa
        ofstream input_compressed(compressed_file, std::ios_base::app);
        ifstream input_unzipped("bcalm.unitigs.shaved.merged.unzipped.gfa");
        string line;
        while (std::getline(input_unzipped, line))
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
        input_unzipped.close();
    }


    //now let's parse the gfa file and decompress it
    kmers = unordered_map<string, pair<string,string>>();
    cout << "Going through the reads again\n";
    go_through_the_reads_again(input_file, "bcalm.unitigs.shaved.merged.unzipped.gfa", context_length, compression, km, kmers);
    cout << "Decompressing\n";
    string decompressed_file = "bcalm.unitigs.shaved.merged.decompressed.gfa";
    expand("bcalm.unitigs.shaved.merged.unzipped.gfa", decompressed_file, km, kmers);

    //test pop_and_shave_homopolymer_errors
    string pop_and_shave_file = "bcalm.unitigs.shaved.merged.decompressed.pop_and_shaved.gfa";
    pop_and_shave_homopolymer_errors(decompressed_file, pop_and_shave_file);
    // exit(0);

    //test compute_eact_cigar
    compute_exact_CIGARs(pop_and_shave_file, output_file, 2000);

    //convert to fasta
    cout << "Convert to fasta\n";
    gfa_to_fasta(output_file, output_file.substr(0, output_file.find_last_of('.')) + ".fasta");

    return 0;
}

