#include "basic_graph_manipulation.h"

#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <algorithm>

using std::cout;
using std::endl;
using std::string;
using std::set;
using std::unordered_map;
using std::cerr;
using std::pair;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::unordered_set;

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
                ss >> depth_string;
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

        if (contig1 == "167--1" || contig2 == "167--1"){
                cout << "found 167--1\n";
            }
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

            int reverse1 = 1; //no reverse complement
            if (std::find(linked[contig1].second.begin(), linked[contig1].second.end(), contig) != linked[contig1].second.end()){
                seq1 = reverse_complement(seq1);
                reverse1 = -1;
            }
            int reverse2 = 1; //no reverse complement
            if (std::find(linked[contig2].second.begin(), linked[contig2].second.end(), contig) != linked[contig2].second.end()){
                seq2 = reverse_complement(seq2);
                reverse2 = -1;
            }

            //now compare the two sequences
            int idx_1 = 0;
            if (reverse1 == -1){
                idx_1 = seq1.size() - 1;
            }
            int idx_2 = 0;
            if (reverse2 == -1){
                idx_2 = seq2.size() - 1;
            }
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
                    idx_1+=reverse1;
                }
                base = seq2[idx_2];
                while (idx_2 < seq2.size() && seq2[idx_2] == base){
                    idx_2+=reverse2;
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

            int reverse1 = 1; //no reverse complement
            if (std::find(linked[contig1].second.begin(), linked[contig1].second.end(), contig) != linked[contig1].second.end()){
                seq1 = reverse_complement(seq1);
                reverse1 = -1;
            }
            int reverse2 = 1; //no reverse complement
            if (std::find(linked[contig2].second.begin(), linked[contig2].second.end(), contig) != linked[contig2].second.end()){
                seq2 = reverse_complement(seq2);
                reverse2 = -1;
            }

            //now compare the two sequences
            int idx_1 = 0;
            if (reverse1 == -1){
                idx_1 = seq1.size() - 1;
            }
            int idx_2 = 0;
            if (reverse2 == -1){
                idx_2 = seq2.size() - 1;
            }
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
                    idx_1+=reverse1;
                }
                base = seq2[idx_2];
                while (idx_2 < seq2.size() && seq2[idx_2] == base){
                    idx_2+=reverse2;
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

void compute_exact_CIGARs(std::string gfa_in, std::string gfa_out, int max_overlap){

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
                cerr << "ERROR: no overlap found between " << name1 << " and " << name2 << "\n";
                exit(1);
            }
            cout << "overlap between " << name1 << " and " << name2 << " is " << overlap << "\n";
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

void create_gaf_from_unitig_graph(std::string unitig_graph, int km, std::string reads_file, std::string gaf_out){
    
    unordered_map<string, string> kmers_to_contigs; //in what contig is the kmer (only unique kmer ofc)
    unordered_map<string, string> reverse_kmers_to_contigs; //what kmers are in the contig

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
            for (int i = 0 ; i <= sequence.size()-km ; i++){
                string kmer = sequence.substr(i, km);
                if (kmers_to_contigs.find(kmer) == kmers_to_contigs.end()){
                    kmers_to_contigs[kmer] = name;
                    reverse_kmers_to_contigs[reverse_complement(kmer)] = name;
                }
                else{
                    kmers_to_contigs[kmer] = "";
                    reverse_kmers_to_contigs[reverse_complement(kmer)] = "";
                }
            }
        }
    }
    input.close();

    //now go through the reads and find the paths
    unordered_map <string, Path> paths;
    ifstream input2(reads_file);
    string name;
    bool nextline = false;
    while (std::getline(input2, line))
    {
        if (line[0] == '@' || line[0] == '>')
        {
            name = line.substr(1, line.size()-1);
            nextline = true;
        }
        else if (nextline){
            Path p;
            //go through the sequence and find the kmers
            for (int i = 0 ; i <= (int)line.size()-km ; i++){
                string kmer = line.substr(i, km);
                if (kmers_to_contigs.find(kmer) != kmers_to_contigs.end() && kmers_to_contigs[kmer] != ""){
                    if (p.contigs.size() > 0 && p.contigs[p.contigs.size()-1] == kmers_to_contigs[kmer] && p.orientations[p.orientations.size()-1] == true){
                        //same contig, do nothing
                    }
                    else{
                        p.contigs.push_back(kmers_to_contigs[kmer]);
                        p.orientations.push_back(true);
                    }

                }
                else if (reverse_kmers_to_contigs.find(kmer) != reverse_kmers_to_contigs.end() && reverse_kmers_to_contigs[kmer] != ""){
                    if (p.contigs.size() > 0 && p.contigs[p.contigs.size()-1] == reverse_kmers_to_contigs[kmer] && p.orientations[p.orientations.size()-1] == false){
                        //same contig, do nothing
                    }
                    else{
                        p.contigs.push_back(reverse_kmers_to_contigs[kmer]);
                        p.orientations.push_back(false);
                    }
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
