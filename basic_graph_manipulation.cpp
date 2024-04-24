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
    system(command.c_str());
}


struct Path{
    vector<string> contigs;
    vector<bool> orientations;
};

/**
 * @brief Create a gaf from unitig graph object and a set of reads. Also compute the coverage of the contigs
 * 
 * @param unitig_graph 
 * @param km 
 * @param reads_file 
 * @param gaf_out
 * @param coverages
 */
void create_gaf_from_unitig_graph(std::string unitig_graph, int km, std::string reads_file, std::string gaf_out, robin_hood::unordered_map<std::string, float>& coverages){
    
    unordered_flat_map<uint64_t, pair<string,int>> kmers_to_contigs; //in what contig is the kmer and at what position (only unique kmer ofc, meant to work with unitig graph)
    unordered_flat_map<string, int> length_of_contigs;

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

            uint64_t hash_foward = -1;
            size_t pos_end = 0;
            long pos_begin = -km;
            // (uint64_t &foward_hash, int k, std::string &seq, size_t &pos )
            while (roll_f(hash_foward, km, sequence, pos_end, pos_begin, false)){
                
                if (pos_begin>=0){
                    kmers_to_contigs[hash_foward] = make_pair(name, pos_end-km);
                }
            }
        }
    }
    input.close();

    //now go through the reads and find the paths
    unordered_flat_map <string, Path> paths;
    ifstream input2(reads_file);
    string name;
    bool nextline = false;
    int nb_reads = 0;
    auto time_now = std::chrono::system_clock::now();
    while (std::getline(input2, line))
    {
        // if (nb_reads%1000==0){
        //     auto time_now2 = std::chrono::system_clock::now();
        //     cout << "aligned " << nb_reads << " on the graph, taking on average " << std::chrono::duration_cast<std::chrono::microseconds>(time_now2 - time_now).count() / (nb_reads+1) << " us per read\r";
        // }
        nb_reads++;
        // cout << "mljqdklmjm " << line << endl;

        if (line[0] == '@' || line[0] == '>')
        {
            name = line.substr(1, line.size()-1);
            nextline = true;
        }
        else if (nextline){
            // if (name.substr(0, name.find_first_of(' ')) != "SRR21295163.5251"){
            //     // cout << "skippple name " << name.substr(0, name.find_first_of(' ')) << "\n";
            //     continue;
            // }
            Path p;
            //go through the sequence and find the kmers
            if (line.size() < km){
                continue;
            }

            uint64_t hash_foward = 0;
            uint64_t hash_reverse = 0;

            int pos_to_look_at = 0;
            size_t pos_end = 0;
            long pos_begin = -km;
            while(roll(hash_foward, hash_reverse, km, line, pos_end, pos_begin, false)){
                if (pos_begin == pos_to_look_at){

                    unsigned long kmer = hash_foward; 

                    if (kmers_to_contigs.find(kmer) != kmers_to_contigs.end()){ //foward kmer
                        if (p.contigs.size() > 0 && p.contigs[p.contigs.size()-1] == kmers_to_contigs[kmer].first && p.orientations[p.orientations.size()-1] == true){
                            //same contig, do nothing
                        }
                        else{
                            string contig = kmers_to_contigs[kmer].first;
                            p.contigs.push_back(contig);
                            p.orientations.push_back(true);
                            if (coverages.find(contig) == coverages.end()){
                                coverages[contig] = 0;
                            }
                            coverages[contig]+= min(1.0, (line.size()- pos_begin) / (double)length_of_contigs[contig]);
                        }
                        //skip the next kmers
                        int length_of_contig_left = length_of_contigs[kmers_to_contigs[kmer].first] - kmers_to_contigs[kmer].second - km;
                        if (length_of_contig_left > 10){ //don't skip too close to the end, you may miss the next contig
                            // pos_to_look_at += (int) length_of_contig_left*0.8; // *0.8 to be sure not to skip the next contig
                            // cout << "dqddf" << endl;
                        }
                        // cout << "found in " << kmers_to_contigs[kmer].first << " at pos " << kmers_to_contigs[kmer].second <<" " << pos_nth << "\n";
                    }
                    else if (kmers_to_contigs.find(hash_reverse) != kmers_to_contigs.end()){ //reverse kmer
                        if (p.contigs.size() > 0 && p.contigs[p.contigs.size()-1] == kmers_to_contigs[hash_reverse].first && p.orientations[p.orientations.size()-1] == false){
                            //same contig, do nothing
                        }
                        else{
                            string contig = kmers_to_contigs[hash_reverse].first;
                            p.contigs.push_back(contig);
                            p.orientations.push_back(false);
                            if (coverages.find(contig) == coverages.end()){
                                coverages[contig] = 0;
                            }
                            coverages[contig]+= min(1.0, (line.size()- pos_begin) / (double)length_of_contigs[contig]);
                        }
                        //skip the next kmers
                        int length_of_contig_left = kmers_to_contigs[hash_reverse].second;
                        if (length_of_contig_left > 10){ //don't skip too close to the end, you may miss the next contig
                            // pos_to_look_at += (int) length_of_contig_left*0.8; // *0.8 to be sure not to skip the next contig
                            // cout << "dqddf" << endl;
                        }
                        // cout << "found in " << kmers_to_contigs[nth.get_reverse_hash()].first << " at pos " << kmers_to_contigs[nth.get_reverse_hash()].second << " " << pos_nth<< "\n";
                    }
                    pos_to_look_at++;
                }
            }

            if (p.contigs.size() > 0){
                paths[name] = p;
            }
        }
    }

    input2.close();

    //now write the gaf file
    ofstream output(gaf_out);
    for (auto p: paths){
        output << p.first << "\t-1\t-1\t-1\t+\t";
        for (int i = 0 ; i < p.second.contigs.size() ; i++){
            if (p.second.orientations[i]){
                output << ">" << p.second.contigs[i];
            }
            else{
                output << "<" << p.second.contigs[i];
            }
        }
        output << "\t\n";
    }

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
        system(convert_command2.c_str());

        //remove tmp files
        string remove_tmp_files = "rm "+path_tmp_folder+"tmp_324*";
        system(remove_tmp_files.c_str());
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
    system(command.c_str());

}

/**
 * @brief Function that takes as input a graph, trim the tips less frequent than abundance_min, remove the branch of bubbles less abundant than abundance_min
 * 
 * @param gfa_in 
 * @param abundance_min
 * @param gfa_out 
 */
void pop_and_shave_graph(string gfa_in, int abundance_min, int min_length, string gfa_out){
    ifstream input(gfa_in);
    //first go through the gfa and find all the places where an end of contig is connected with two links
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

            char or1 = (orientation1 == "+" ? 1 : 0);
            char or2 = (orientation2 == "+" ? 0 : 1);
            if (orientation1 == "+"){
                linked[name1].second.push_back(make_pair(name2, or2));
            }
            else{
                linked[name1].first.push_back(make_pair(name2, or2));
            }

            if (orientation2 == "+"){
                linked[name2].first.push_back(make_pair(name1, or1));
            }
            else{
                linked[name2].second.push_back(make_pair(name1, or1));
            }
        }
        pos += line.size() + 1;
    }
    input.close();
    input.open(gfa_in);

    unordered_set<string> to_keep; //kept contigs are the one with a coverage above abundance_min and their necessary neighbors for the contiguity

    //iterative cleaning of the graph
    int number_of_edits = 1;
    while (number_of_edits > 0){
        number_of_edits = 0;

        for (auto c: linked){

            string contig = c.first;

            if (coverage[contig] >= abundance_min || length_of_contigs[contig] >= min_length){
                to_keep.insert(contig);
            }

            if (to_keep.find(contig) != to_keep.end()){
                //make sure the contig has at least one neighbor left and right (if not, take the one with the highest coverage)
                float best_coverage = 0;
                string best_contig = "";
                for (auto l: linked[contig].second){
                    if (coverage[l.first] > best_coverage){
                        best_coverage = coverage[l.first];
                        best_contig = l.first;
                    }
                }
                if (best_contig != ""){
                    if (to_keep.find(best_contig) == to_keep.end()){
                        to_keep.insert(best_contig);
                        number_of_edits++;
                    }
                }
            }
        }

        //update the links by removing all references to the contigs to remove

        // cout << "Finished one iteration of pop and shave, added " << number_of_edits << " contigs to keep\n";

    }

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

            if (to_keep.find(name1) != to_keep.end() && to_keep.find(name2) != to_keep.end()){
                out << line << "\n";
            }
        }
    }
    out.close();

}