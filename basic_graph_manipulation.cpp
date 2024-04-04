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


void pop_and_shave_homopolymer_errors(string gfa_in, string gfa_out){
    ifstream input(gfa_in);
    //first go through the gfa and find all the places where an end of contig is connected with two links
    unordered_map<string, pair<vector<string>, vector<string>>> linked;
    unordered_map<string, long int> pos_of_contig_seq_in_file;
    unordered_map<string, float> coverage;

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

            if (linked.find(name1) == linked.end()){
                linked[name1] = {vector<string>(0), vector<string>(0)};
            }
            if (linked.find(name2) == linked.end()){
                linked[name2] = {vector<string>(0), vector<string>(0)};
            }
            if (orientation1 == "+"){
                linked[name1].second.push_back(name2);
            }
            else{
                linked[name1].first.push_back(name2);
            }

            if (orientation2 == "+"){
                linked[name2].first.push_back(name1);
            }
            else{
                linked[name2].second.push_back(name1);
            }
        }
        pos += line.size() + 1;
    }
    input.close();
    input.open(gfa_in);

    unordered_set<string> to_remove;

    //now find all the places where there is a node
    for (auto c: linked){
        string contig = c.first;
        if (c.second.first.size() == 2){

            string contig1 = c.second.first[0];
            string contig2 = c.second.first[1];

            // if (contig1 == "167--1" || contig2 == "167--1"){
            //     cout << "found 167--1\n";
            // }
            bool bubble = false; //either a bubble or a dead end
            //check if there is a dead end or a bubble
            if (linked[contig1].first.size() == 0 || linked[contig1].second.size() == 0 || linked[contig2].first.size() == 0 || linked[contig2].second.size() == 0){
                //dead end, we love it
            }
            else if (linked[contig1].first.size() == 1 && linked[contig1].second.size() == 1 && linked[contig2].first.size() == 1 && linked[contig2].second.size() == 1
                && ((linked[contig1].first[0] == linked[contig2].first[0] && linked[contig1].second[0] == linked[contig2].second[0]) ||
                    linked[contig1].second[0] == linked[contig2].first[0] && linked[contig1].first[0] == linked[contig2].second[0]) ){
                //bubble, we love it
                bubble = true;
            }
            else{
                //something else, we don't love it, keep the contigs
                continue;
            }

            //compare the two contigs and see if they are identical homopolymer-wise or not
            input.seekg(pos_of_contig_seq_in_file[contig1]);
            std::getline(input, line);
            string dont_care, seq1, seq2;
            std::stringstream ss(line);
            ss >> dont_care >> dont_care >> seq1;
            input.seekg(pos_of_contig_seq_in_file[contig2]);
            std::getline(input, line);

            ss = std::stringstream(line);
            ss >> dont_care >> dont_care >> seq2;

            if (std::find(linked[contig1].second.begin(), linked[contig1].second.end(), contig) != linked[contig1].second.end()){
                seq1 = reverse_complement(seq1);
            }
            int reverse2 = 1; //no reverse complement
            if (std::find(linked[contig2].second.begin(), linked[contig2].second.end(), contig) != linked[contig2].second.end()){
                seq2 = reverse_complement(seq2);
            }

            //now compare the two sequences
            int idx_1 = 0;
            int idx_2 = 0;
            bool identical = true;
            while (idx_1 < seq1.size() && idx_2 < seq2.size()){
                if (seq1[idx_1] != seq2[idx_2]){
                    identical = false;
                    // if (contig1 == "167--1" || contig2 == "167--1"){
                    //     cout << "not identical at " << idx_1 << " and " << idx_2 << "\n";
                    //     cout << "seq1 " << seq1.substr(0, idx_1) << " " << seq1.substr(idx_1+1, std::min((int)seq1.size()-idx_1-1, 100)) << "\n";
                    //     cout << "seq2 " << seq2.substr(0, idx_2) << " " << seq2.substr(idx_2+1, seq2.size()-idx_2-1) << "\n";
                    // }
                    break;
                }
                  
                char base = seq1[idx_1];
                while (idx_1 < seq1.size() && idx_1 >= 0 && seq1[idx_1] == base){
                    idx_1++;
                }
                base = seq2[idx_2];
                while (idx_2 < seq2.size() && idx_2 >= 0 && seq2[idx_2] == base){
                    idx_2++;
                }
            }

            if (identical){
                cout << "contig " << contig1 << " and " << contig2 << " are identical\n";
                int coverage1 = coverage[contig1];
                int coverage2 = coverage[contig2];
                int length1 = seq1.size();
                int length2 = seq2.size();

                if (coverage1 > coverage2 && (bubble || length1 > length2)){
                    // //keep contig1
                    // cout << "keeping contig " << contig1 << "\n";
                    // //remove contig2
                    // cout << "removing contig " << contig2 << "\n";
                    to_remove.insert(contig2);
                }
                else if (coverage2 >= coverage1 && (bubble || length2 > length1)){
                    //keep contig2
                    // cout << "keeping contig " << contig2 << "\n";
                    // //remove contig1
                    // cout << "removing contig " << contig1 << "\n";
                    to_remove.insert(contig1);
                }
                else{
                    cout << "ERROR ambiguous call here between " << contig1 << " and " << contig2 << "\n"; //the dead end 
                    // exit(1);
                }
            }
        }
        if (c.second.second.size() == 2){

            string contig1 = c.second.second[0];
            string contig2 = c.second.second[1];

            // if (contig1 == "167--1" || contig2 == "167--1"){
            //     cout << "found 167--1\n";
            // }
            bool bubble = false; //either a bubble or a dead end
            //check if there is a dead end or a bubble
            if (linked[contig1].first.size() == 0 || linked[contig1].second.size() == 0 || linked[contig2].first.size() == 0 || linked[contig2].second.size() == 0){
                //dead end, we love it
            }
            else if (linked[contig1].first.size() == 1 && linked[contig1].second.size() == 1 && linked[contig2].first.size() == 1 && linked[contig2].second.size() == 1
                && ((linked[contig1].first[0] == linked[contig2].first[0] && linked[contig1].second[0] == linked[contig2].second[0]) ||
                    linked[contig1].second[0] == linked[contig2].first[0] && linked[contig1].first[0] == linked[contig2].second[0]) ){
                //bubble, we love it
                bubble = true;
            }
            else{
                //something else, we don't love it, keep the contigs
                continue;
            }

            //compare the two contigs and see if they are identical homopolymer-wise or not
            input.seekg(pos_of_contig_seq_in_file[contig1]);
            std::getline(input, line);
            string dont_care, seq1, seq2;
            std::stringstream ss(line);
            ss >> dont_care >> dont_care >> seq1;
            input.seekg(pos_of_contig_seq_in_file[contig2]);
            std::getline(input, line);

            ss = std::stringstream(line);
            ss >> dont_care >> dont_care >> seq2;

            if (std::find(linked[contig1].second.begin(), linked[contig1].second.end(), contig) != linked[contig1].second.end()){
                seq1 = reverse_complement(seq1);
            }
            if (std::find(linked[contig2].second.begin(), linked[contig2].second.end(), contig) != linked[contig2].second.end()){
                seq2 = reverse_complement(seq2);
            }

            //now compare the two sequences
            int idx_1 = 0;
            int idx_2 = 0;
            bool identical = true;
            while (idx_1 < seq1.size() && idx_2 < seq2.size()){
                if (seq1[idx_1] != seq2[idx_2]){
                    identical = false;
                    // if (contig1 == "167--1" || contig2 == "167--1"){
                    //     cout << "not identical at " << idx_1 << " and " << idx_2 << "\n";
                    //     cout << "seq1 " << seq1.substr(0, idx_1) << " " << seq1.substr(idx_1+1, std::min((int)seq1.size()-idx_1-1, 100)) << "\n";
                    //     cout << "seq2 " << seq2.substr(0, idx_2) << " " << seq2.substr(idx_2+1, seq2.size()-idx_2-1) << "\n";
                    // }
                    break;
                }
                  
                char base = seq1[idx_1];
                while (idx_1 < seq1.size() && seq1[idx_1] == base){
                    idx_1++;
                }
                base = seq2[idx_2];
                while (idx_2 < seq2.size() && seq2[idx_2] == base){
                    idx_2++;
                }
            }

            if (identical){
                cout << "contig " << contig1 << " and " << contig2 << " are identical\n";
                int coverage1 = coverage[contig1];
                int coverage2 = coverage[contig2];
                int length1 = seq1.size();
                int length2 = seq2.size();

                if (coverage1 > coverage2){
                    // //keep contig1
                    // cout << "keeping contig " << contig1 << "\n";
                    // //remove contig2
                    // cout << "removing contig " << contig2 << "\n";
                    to_remove.insert(contig2);
                }
                else if (coverage2 >= coverage1){
                    //keep contig2
                    // cout << "keeping contig " << contig2 << "\n";
                    // //remove contig1
                    // cout << "removing contig " << contig1 << "\n";
                    to_remove.insert(contig1);
                }
                else{
                    cout << "ambiguous call here between " << contig1 << " and " << contig2 << "\n"; //the dead end 
                    exit(1);
                }
            }
        }
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
            if (to_remove.find(name) == to_remove.end()){
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

            if (to_remove.find(name1) == to_remove.end() && to_remove.find(name2) == to_remove.end()){
                out << line << "\n";
            }
        }
    }
    out.close();

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
        if (nb_reads%1000==0){
            auto time_now2 = std::chrono::system_clock::now();
            cout << "aligned " << nb_reads << " on the graph, taking on average " << std::chrono::duration_cast<std::chrono::microseconds>(time_now2 - time_now).count() / (nb_reads+1) << " us per read\r";
        }
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


void merge_adjacent_contigs_BCALM(std::string gfa_in, std::string gfa_out, int k, std::string path_to_bcalm, std::string path_convertToGFA){
        
        //convert gfa_in to fasta
        cout << "Convert to fasta bcalm.unitigs.shaved.gfa\n";
        string tmp_fasta = "tmp_324.fasta";
        gfa_to_fasta(gfa_in, tmp_fasta);

        //to merge, simply make a unitig graph from bcalm.unitigs.shaved.gfa and then convert it to gfa
        cout << "Creating shaved unitig graph\n";
        string command_unitig_graph = path_to_bcalm + " -in " + tmp_fasta + " -kmer-size "+std::to_string(k)+" -abundance-min 1 -out tmp_324 > bcalm.log 2>&1";
        auto unitig_graph_ok = system(command_unitig_graph.c_str());
        cout << "launching unitig graph\n" << command_unitig_graph << endl;
        if (unitig_graph_ok != 0){
            cerr << "ERROR: unitig graph failed in merge_adjacent_contigs_BCALM\n";
            cout << command_unitig_graph << endl;
            exit(1);
        }

        //convert to gfa
        cout << "Launching convertToGFA\n";
        string convert_command2 = path_convertToGFA + " tmp_324.unitigs.fa " + gfa_out + " " + std::to_string(k);
        system(convert_command2.c_str());

        //remove tmp files
        string remove_tmp_files = "rm tmp_324*";
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
 * @brief When two contigs are identical homopolymer-wise, keep the one with the highest coverage, and add the two coverages
 * 
 * @param gfa_in 
 * @param gfa_out 
 */
void remove_homopolymer_errors(std::string gfa_in, std::string gfa_out, std::string path_to_minimap2){

    //first create a temporary homopolymer-compressed version of gfa_in. Also index the coverages of the contigs and their number of neighbors
    string tmp_gfa = gfa_in + ".hpc.tmp.fa";
    unordered_map<string, float> coverages;
    unordered_map<string, pair<set<string>,set<string>>> neighbors;
    ifstream input(gfa_in);
    ofstream out(tmp_gfa);
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
            string hpc_sequence = "";
            char base = sequence[0];
            for (int i = 1 ; i < sequence.size() ; i++){
                if (sequence[i] != base){
                    hpc_sequence += base;
                    base = sequence[i];
                }
            }
            hpc_sequence += base;
            out << ">" << name << "\n" << hpc_sequence << "\n";

            //now index the coverage
            string depth_string = "  ";
            while (depth_string.size()>= 2 && (depth_string.substr(0,5) != "DP:f:")){
                string nds;
                ss >> nds;
                depth_string = nds;
            }
            if (depth_string.substr(0,5) == "DP:f:"){
                float depth = std::stof(depth_string.substr(5, depth_string.size()-5));
                coverages[name] = depth;
            }
            else{
                coverages[name] = 0;
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

            if (neighbors.find(name1) == neighbors.end()){
                neighbors[name1] = {set<string>(), set<string>()};
            }
            if (neighbors.find(name2) == neighbors.end()){
                neighbors[name2] = {set<string>(), set<string>()};
            }
            if (orientation1 == "+"){
                neighbors[name1].second.emplace(name2);
            }
            else{
                neighbors[name1].first.emplace(name2);
            }
            if (orientation2 == "+"){
                neighbors[name2].first.emplace(name1);
            }
            else{
                neighbors[name2].second.emplace(name1);
            }
        }
    }
    out.close();
    input.close();

    //now align the contigs to each other with minimap2 asking for 100% identity
    string tmp_paf = gfa_in + ".hpc.tmp.paf";
    string command = path_to_minimap2 + " -x asm5 -k 25 -X -c -B 10 -O 10 -E 10 " + tmp_gfa + " " + tmp_gfa + " > " + tmp_paf;
    auto minimap_res = system(command.c_str());
    if (minimap_res != 0){
        cerr << "ERROR: minimap2 failed in remove_homopolymer_errors\n";
        cout << command << endl;
        exit(1);
    }

    //now go through the paf file and find the contigs that are identical homopolymer-wise
    ifstream input2(tmp_paf);
    std::set<string> toRemove;
    while (std::getline(input2, line))
    {
        std::stringstream ss(line);
        string name1, name2, dont_care;
        int length1, pos1_1, pos1_2, length2, pos2_1, pos2_2;
        ss >> name1 >> length1 >> pos1_1 >> pos1_2 >> dont_care >> name2 >> length2 >> pos2_1 >> pos2_2;
        if (name1 == name2 || ((pos1_2-pos1_1 != length1) && (pos2_2-pos2_1 != length2))){
            continue;
        }
        //now look for the de:f: tag and check that it is de:f:0
        string de_tag;
        while (ss >> de_tag){
            if (de_tag.substr(0, 5) == "de:f:"){
                if (de_tag == "de:f:0"){

                    float coverage1 = coverages[name1];
                    float coverage2 = coverages[name2];

                    //now check the configurations

                    //1st configuration: one contig is fully contained in the other and has no neighbors or is a dead end, delete it
                    if (pos1_1 == 0 && pos1_2 == length1 && (neighbors[name1].first.size() == 0 || neighbors[name1].second.size() == 0)){
                        toRemove.insert(name1);
                    }
                    else if (pos2_1 == 0 && pos2_2 == length2 && (neighbors[name2].first.size() == 0 || neighbors[name2].second.size() == 0)){
                        toRemove.insert(name2);
                    }

                    //2nd configuration: the two contigs are identical and form a bubble, keep the one with the highest coverage
                    if (neighbors[name1].first == neighbors[name2].first && neighbors[name1].second == neighbors[name2].second){
                        if (coverage1 > coverage2){
                            toRemove.insert(name2);
                        }
                        else{
                            toRemove.insert(name1);
                        }
                    }
                    else if (neighbors[name1].first == neighbors[name2].second && neighbors[name1].second == neighbors[name2].first){
                        if (coverage1 > coverage2){
                            toRemove.insert(name2);
                        }
                        else{
                            toRemove.insert(name1);
                        }
                    }
                }
                break;
            }
        }
    }

    //now go through the gfa file and remove the contigs
    input.open(gfa_in);
    out.open(gfa_out);

    cout << "Remoovinng..." << endl;
    for (auto c: toRemove){
        cout << c << endl;
    }

    while (std::getline(input, line))
    {
        if (line[0] == 'S')
        {
            string name;
            string dont_care;
            string sequence;
            std::stringstream ss(line);
            ss >> dont_care >> name >> sequence;
            if (toRemove.find(name) == toRemove.end()){
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

            if (toRemove.find(name1) == toRemove.end() && toRemove.find(name2) == toRemove.end()){
                out << line << "\n";
            }
        }
    }
    out.close();

    //remove tmp files
    string remove_tmp_files = "rm " + tmp_gfa + " " + tmp_paf;
    system(remove_tmp_files.c_str());
}


