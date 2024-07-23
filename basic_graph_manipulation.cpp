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
    ifstream input2(reads_file);
    string name;
    bool nextline = false;
    int nb_reads = 0;
    auto time_now = std::chrono::system_clock::now();
    ofstream output(gaf_out);
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
                            pos_to_look_at += (int) length_of_contig_left*0.8; // *0.8 to be sure not to skip the next contig
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
                            pos_to_look_at += (int) length_of_contig_left*0.8; // *0.8 to be sure not to skip the next contig
                            // cout << "dqddf" << endl;
                        }
                        // cout << "found in " << kmers_to_contigs[nth.get_reverse_hash()].first << " at pos " << kmers_to_contigs[nth.get_reverse_hash()].second << " " << pos_nth<< "\n";
                    }
                    pos_to_look_at++;
                }
            }

            if (p.contigs.size() > 0){
                output << name << "\t-1\t-1\t-1\t+\t";
                for (int i = 0 ; i < p.contigs.size() ; i++){
                    if (p.orientations[i]){
                        output << ">" << p.contigs[i];
                    }
                    else{
                        output << "<" << p.contigs[i];
                    }
                }
                output << "\t\n";}
        }
    }

    input2.close();
    output.close();

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
    system(command.c_str());

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
 * @return int 
 */
int explore_neighborhood(string contig, int endOfContig, unordered_map<string, pair<vector<pair<string, char>>, vector<pair<string,char>>>> &linked, unordered_map<string, float>& coverage, unordered_map<string, int>& length_of_contigs, int k, int length_left, float original_coverage, unordered_map<string, pair<bool, bool>> &not_overcovered){
    
    if (length_left <= 0){
        return 2;
    }

    bool overcovered_in_neighborhood = false;
    if (endOfContig == 1){ //arrived by the right, go through to the left
        if (linked[contig].first.size() == 0){
            return 0;
        }
        for (auto l: linked[contig].first){
            if (coverage[l.first] > 20*original_coverage){
                overcovered_in_neighborhood = true;
            }
            else if (l.first != contig) { //the condition is so we don't end up in a loop
                if (not_overcovered.find(l.first) != not_overcovered.end() && not_overcovered[l.first].first && l.second == 1 && coverage[l.first] < 2*coverage[contig] && coverage[l.first]*2 > original_coverage){
                    return 2;
                }
                else if (not_overcovered.find(l.first) != not_overcovered.end() && not_overcovered[l.first].second && l.second == 0 && coverage[l.first] < 2*coverage[contig] && coverage[l.first]*2 > original_coverage){
                    return 2;
                }

                int res = explore_neighborhood(l.first, l.second, linked, coverage, length_of_contigs, k, length_left - length_of_contigs[l.first] + k - 1, original_coverage, not_overcovered);
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
            if (coverage[l.first] > 20*original_coverage){
                overcovered_in_neighborhood = true;
            }
            else {

                if (not_overcovered.find(l.first) != not_overcovered.end() && not_overcovered[l.first].first && l.second == 1 && coverage[l.first] < 2*original_coverage && coverage[l.first]*2 > original_coverage){
                    return 2;
                }
                else if (not_overcovered.find(l.first) != not_overcovered.end() && not_overcovered[l.first].second && l.second == 0 && coverage[l.first] < 2*original_coverage && coverage[l.first]*2 > original_coverage){
                    return 2;
                }

                int res = explore_neighborhood(l.first, l.second, linked, coverage, length_of_contigs, k, length_left - length_of_contigs[l.first] + k - 1, original_coverage, not_overcovered);
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
 * @param abundance_min
 * @param gfa_out 
 */
void pop_and_shave_graph(string gfa_in, int abundance_min, int min_length, bool contiguity, int k, string gfa_out){
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
    int number_of_edits = 1;
    unordered_map<string, pair<bool, bool>> not_overcovered; //iteratively mark the contigs that are not overcovered left or right
    not_overcovered.reserve(linked.size());
    while (number_of_edits > 0){
        number_of_edits = 0;
        int count = 0;

        for (auto c: linked){

            string contig = c.first;
            // cout << "contig " << contig << " djuqst sdkdqsm \n";
            // if (contig != "254117"){
            //     continue;
            // }

            if (coverage[contig] >= abundance_min || length_of_contigs[contig] >= min_length){
                //this contig is solid EXCEPT if it is a tip that has an abundance less than 10% of the coverage of the contig it is connected to
                if (linked[contig].first.size() == 0 && linked[contig].second.size() == 1){
                    //this is a tip
                    if (coverage[contig] >= 0.1*coverage[linked[contig].second[0].first]){
                        to_keep.insert(contig);
                    }
                }
                else if (linked[contig].first.size() == 1 && linked[contig].second.size() == 0){
                    //this is a tip
                    if (coverage[contig] >= 0.1*coverage[linked[contig].first[0].first]){
                        to_keep.insert(contig);
                    }
                }
                else if (contiguity) {//if contiguity mode is on, do not take contigs that have in all their neighborhood contigs with 20x their coverage (corresponds to poorly covered side of a bubble)
                    
                    int size_of_neighborhood = 5*k;
                    // cout << "launching..\n";
                    int overcovered_right = 2;
                    if (not_overcovered.find(contig) == not_overcovered.end() || !not_overcovered[contig].second){
                        overcovered_right = explore_neighborhood(contig, 0, linked, coverage, length_of_contigs, k, size_of_neighborhood, coverage[contig], not_overcovered);
                    }
                    if (overcovered_right == 2){
                        if (not_overcovered.find(contig) != not_overcovered.end()){
                            not_overcovered[contig].first = false;
                            not_overcovered[contig].second = false;
                        }
                        not_overcovered[contig].second = true;
                    }
                    int overcovered_left = 2;
                    if (not_overcovered.find(contig) == not_overcovered.end() || !not_overcovered[contig].first){
                        overcovered_left = explore_neighborhood(contig, 1, linked, coverage, length_of_contigs, k, size_of_neighborhood, coverage[contig], not_overcovered);
                    }
                    if (overcovered_left == 2){
                        if (not_overcovered.find(contig) != not_overcovered.end()){
                            not_overcovered[contig].first = false;
                            not_overcovered[contig].second = false;
                        }
                        not_overcovered[contig].first = true;
                    }
                    // cout << contig << " " << overcovered_left << " " << overcovered_right << "\n";

                    // if (contig == "280245"){
                    //     cout << "contig " << contig << " has overcovered neighborhood " << overcovered_left << " " << overcovered_right << " " << coverage[contig] << "\n";
                    //     exit(1);
                    // }

                    if ((overcovered_left == 2 || overcovered_right == 2) || (overcovered_left == 0 && overcovered_right == 0)){
                        to_keep.insert(contig);
                    }
                    // else{
                    //     cout << "contig " << contig << " has overcovered neighborhood " << count << " " << linked.size() << "\n";
                    //     // exit(1);
                    // }
                }
                else{
                    to_keep.insert(contig);
                }
                count ++;
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

/**
 * @brief Recursive function that lists all the paths from a given end of a contig up to a certain distance
 * 
 * @param contig_name 
 * @param contig_end 
 * @param length_to_explore 
 * @param linked 
 * @param length_of_contigs 
 * @return vector<vector<pair<int,int>>> 
 */
vector<vector<pair<string,char>>> list_all_paths_from_this_end(string contig_name, int contig_end, int length_to_explore, unordered_map<string, pair<vector<pair<string, char>>, vector<pair<string,char>>>> &linked, unordered_map<string, int> &length_of_contigs){
    if (length_to_explore <= 0){
        return {{{contig_name, contig_end}}};
    }

    vector<vector<pair<string,char>>> paths;
    if (contig_end == 0){
        for (auto l: linked[contig_name].first){
            vector<vector<pair<string,char>>> paths_from_this_neighbor = list_all_paths_from_this_end(l.first, 1-l.second, length_to_explore - length_of_contigs[l.first], linked, length_of_contigs);
            for (vector<pair<string,char>> p: paths_from_this_neighbor){
                p.push_back({contig_name, 0});
                paths.push_back(p);
            }
        }
    }
    else if (contig_end == 1){
        for (auto l: linked[contig_name].second){
            vector<vector<pair<string,char>>> paths_from_this_neighbor = list_all_paths_from_this_end(l.first, 1-l.second, length_to_explore - length_of_contigs[l.first], linked, length_of_contigs);
            for (vector<pair<string,char>> p: paths_from_this_neighbor){
                p.push_back({contig_name, 1});
                paths.push_back(p);
            }
        }
    }
    return paths;

}

/**
 * @brief In small bubbles, take only one side and discard the other
 * 
 * @param gfa_in 
 * @param length_of_longest_read only pop bubbles that are between two contigs that are at least this length and at most this length/2 apart
 * @param gfa_out 
 */
void pop_bubbles(std::string gfa_in, int length_of_longest_read, std::string gfa_out){

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

    //now pop the bubbles

    std::set<pair<pair<string,char>,pair<string,char>>> links_to_delete;
    //iterate through all contigs: if it is longer than the longest read, see if there is a bubble left or right
    for (auto c: linked){
        string contig_name = c.first;
        // if (contig_name != "227446_461431_221297_221306_526669_119271_119327_121126_119302_121134_119288_534150_220446_431651_667036_319694_319670_319708_319687_319702_675626_323769_323791_323775_323786_662373_217425_217437_217444_456904_271066_271077_271058_737991_57172_57169_57167_57156_668852_513523_602974_29226_521082_7257_683000_684155_265009_264999_265007_265029_265012_265015_265024_713618_717586_359269_359272_468784_120994_121018_509741_70028_70024_506198_9398_9362_9378_9386_513877_748823_221229_454640_322458_322471_322451_717094_331085_720583_732326_7940_695958_17340_593441_18873_18854_758111_597621_14868_477723_275257_275271_275277_275298_440876_131606_600592_540747_223989_631563_493122_541558_73092_73131_73119_73088_73113_465651_363871_363863_363842_363851_363856_363867_363833_592416_366363_366411_664304_318446_733766_220362_220374_220309_220304_220326_220312_220370_504345_172596_441602_370138_370101_370122_370095_370171_370107_370162_370166_651806_619350_130455_566710_327348_327305_327339_327300_327328_564721_385883_478348_230566_230570_230550_583463_276702_276706_276697_276717_276676_610393_281971_517558_696341_377152_377147_377136_377128_377098_734543_383647_383629_383621_531366_35452_35440_669605_754666_76303_76294_76312_76320_76285_76326_631931_372254_372244_372222_372214_372272_435517_136809_728874_230999_231035_231024_231013_231062_231050_655472_579972_281741_281736_281747_281754_281731_477849_188969_188982_188991_532223_476817_29567_29576_29590_29571_508720_720269_137196_137204_462054_29440_29417_29427_29420_29447_29450_764150_73884_73923_73905_73919_73915_73933_73872_572602_76264_76235_76270_76228_667831_279247_279243_279295_279269_279239_438071_335765_335757_335725_335736_335733_335809_335788_628759_84653_84668_424088_333586_333599_333570_333577_602860_725254_318831_318839_701871_16041_16046_16060_16053_442725_278197_278176_278188_701947_222540_222488_222519_222509_222496_441072_221148_221156_221144_221129_757342_272994_273030_273013_273017_722989_71096_71093_71103_71155_71184_71132_71149_425754_365079_365123_365101_365130_365095_748628_548579_73062_73066_73059_686648_274821_274849_274831_274836_274842_490435_20768_20792_20762_20797_20773_756835_368846_368881_368863_368872_735652_489803_618745_361824_361836_361846_361821_455441_223074_223077_223083_732303_20221_20230_714596_225020_224961_224974_224968_224996_225016_224979_224982_461507_17625_17640_741230_376552_376585_750556_68008_67964_67956_67940_67978_67946_684532_119402_119431_119380_119376_119390_435887_70109_70104_70092_731081_57544_57513_57522_57502_57483_57506_57570_502383_212289_212285_624785_378009_378000_377982_696798_231856_231862_436907_223051_223038_745270_238668_238685_238678_482502_28768_28805_28772_28720_28812_28754_28777_493369_755627_329322_722862_181054_181059_763375_291101_291106_736588_326579_326599_326607_658595_224903_224913_421945_488516_533658_68821_68816_68807_68800_68813_68825_456907_365764_560945_500423_325135_325148_325142_542926_284082_681254_170210_170220_170216_170240_170230_545848_20629_571376_122956_122947_526580_272263_272272_272282_586220_374575_679761_489325_30068_546086_34027_34047_474759_273544_273555_273524_273540_273534_273501_579291_126589_457357_185603_593370_397086_125733_125744_429141_274680_274698_274671_274667_274708_274688_633804_168189_168193_612609_275400_275410_551833_270432_270476_270446_270452_270421_515973_426991_7099_7094_7108_7112_7088_7117_587750_77888_77885_77895_754936_19849_674012_18836_18856_18832_581156_466061_175183_175194_175177_175207_175191_175218_764141_29436_29452_29455_453621_188055_188051_395361_86695_86689_86678_663524_373999_373974_673568_376502_376521_572275_681343_28187_508430_273551_273584_273558_273576_558742_598860_181932_181940_519999_468283_607000_553407_329319_329329_511157_173138_173147_547306_752437_17505_17518_17525_627125_176988_592646_123835_123849_123896_123857_123872_123870_123826_769246_131653_131658_131631_131643_414713_227062_608050_30182_30197_30187_407284_364939_364929_596490_224281_224240_224242_224250_224232_224273_224284_224282_224255_224244_563254_632376_666308_170217_170235_170268_170256_170228_741485_74689_74677_494035_75058_507654_323044_323064_494491_684746_234023_234012_536820_18585_18578_585909_217885_217892_217875_619427_279697_734900_322745_322735_745484_312738_312743_640133_15561_15565_15571_425317_274295_274303_274337_274309_274321_274330_721231_274693_274666_274659_274642_274678_274702_685879_324073_324033_324076_324077_324042_324036_324063_324067_449397_225227_608546_165917_165933_165968_165958_165941_562389_19180_647994_326711_482643_164606_164622_164641_164609_164614_164625_164653_475784_265841_265835_474688_166706_166714_166696_166728_166723_166734_166689_394626_567402_275921_275881_275940_275869_275948_275890_275916_706488_173795_173830_173851_173859_173866_740000_461047_179528_179502_735702_369771_367605_367627_369774_369790_529573_216419_216446_216407_216394_216381_216401_216440_216387_414672_379474_729748_126945_126885_126890_126903_126932_754220_7785_7779_7794_7801_432720_289537_590225_237777_237780_237790_700582_276668_686089_326185_326227_326206_326198_326191_546359_643890_68553_677090_376449_463368_690876_222464_222435_222430_222441_222453_435504_21457_21461_675658_368131_368126_654461_280250_280258_607729_125996_126031_126001_126047_126018_691961_275115_275111_402263_607492_595767_278017_760160_173945_173882_173950_173919_173955_173943_173930_543137_314443_314472_312917_314454_312906_314450_314439_665254_213584_213587_213576_213569_581696_329587_603553_322510_322498_322495_322493_538680_13696_13711_13754_13748_13724_13732_13686_712587_613082_14921_14929_696189_311103_311089_311100_311095_670422_364269_526566_209460_209478_209469_209482_209464_514092_54287_54281_440701_520184_14289_14295_12398_14301_12394_14309_14342_747600_548363_60375_60354_60360_60333_60329_60344_422147_659374_374595_374619_374633_374601_374630_469477_227518_227463_227488_227505_744055_168130_168143_432376_117684_640911_600671_10881_10872_10931_10898_10935_10867_10911_481508_322168_322178_453268_713636_269755_269750_269751_639062_472394_171801_171831_171813_171806_171821_171837_614746_323459_323468_435396_364866_364846_628083_436525_268118_268076_268131_268140_268110_268093_423630_725191_626017_323628_323620_323657_323649_323643_323677_323670_408917_114410_114400_114419_114396_114393_114415_398372_691567_71563_718769_328873_328870_627987_185495_185487_588960_7110_7115_7129_7124_7134_677837_448419_763406_23109_23117_717454_94245_94259_602154_213679_213706_213669_213697_213701_213692_761053_176164_176152_176156_176143_176182_432223_621564_168645_168614_168618_168641_168650_533947_494580_663388_311330_311312_311288_311300_311293_311337_311308_311342_721952_216358_739574_272532_272556_488442_59213_59199_59182_59192_491754_654557_363514_363518_363561_363528_363546_363564_363551_363489_659145_275317_275319_659042_367568_367558_618004_164980_164923_164946_164963_164930_164942_164902_481907_21520_21491_736523_734405_172652_172677_172662_471888_322367_322375_322361_552844_225579_225613_225619_486296_325455_325464_325460_411736_686236_708311_168181_510323_456838_14131_14166_14119_14208_14222_14194_14214_445374_681922_402226_122035_122052_669045_60811_60792_60816_60807_60801_704235_177788_433953_364890_364876_364871_448341_13727_748539_321641_321638_321600_321621_321635_661186_123758_123745_123751_123792_123778_726134_127649_667224_180398_739202_16453_16410_16421_16445_16449_16425_735043_317719_317704_317712_655624_369838_369810_369779_369802_369775_369834_766994_634008_759936_320028_320016_320024_671823_13836_13820_13852_13825_13862_13828_13858_481866_232796_232771_232790_504008_484558_125139_125163_125072_125131_125109_125157_454488_231599_231575_231594_726240_333121_333104_333115_333110_333088_333060_436624_10027_741617_223282_223268_223260_713634_753540_721308_129153_129170_129158_129202_129219_769918_285922_685584_69993_70002_69960_69972_69988_69956_410228_371765_371790_371749_371781_642986_180168_180181_180160_180172_645891_526586_127636_127613_127598_127628_459065_619015_231672_231682_231689_231665_531477_226490_226503_661558_126509_126542_126484_126489_126502_126525_126564_126558_457762_173054_173030_173047_173043_173066_611551_332229_332242_538443_171689_171695_643506_223023_223017_533155_637603_322358_322418_322378_322402_322396_448721_524111_284023_284036_284017_284030_617364_517751_365578_730939_171041_171035_171094_415620_422003_237585_683175_375345_375338_543346_132354_732346_328170_328163_328157_328179_448639_662628_25591_25625_25612_25584_577917_319542_319552_319540_319575_319556_319564_743228_617369_187498_187489_187480_187502_187471_187510_514879_228112_228098_228104_407402_76197_76216_76207_76234_548755_25420_25414_25452_25425_690140_382871_382866_463142_338087_338081_338055_338078_338073_544307_335862_562469_187768_187783_484934_693389_451375_376080_376107_376087_495996_132960_457656_368407_368423_368410_368416_368414_432024_226772_442422_234681_599124_339625_339611_572933_518283_228942_228927_412562_22780_22798_22771_414325_35855_35847_35863_699561_71976_71986_71981_72003_71964_71989_430140_427562_127053_127035_448706_120497_458021_175399_599879_464485_15519_15543_15532_15575_15529_15500_15522_15588_733616_642717_365522_365543_365550_529482_24151_24179_24159_482504_73573_756227_328749_621313_174871_174844_174851_759700_284067_284069_612180_679201_123646_123607_123594_123618_123640_123566_123631_123589_123550_516295_744632_766268_120827_120822_120801_120816_696430_221881_221871_221907_221901_221894_221933_221974_441680_476761_22805_22846_22816_22825_589576_171020_171005_170990_170982_526998_375065_374997_375042_375036_375003_375069_405416_171424_171418_638531_510190_323060_323025_323043_323038_323019_323049_693100_726082_14586_405575_165211_165202_165183_165174_165225_671198_127280_499344_543766_722158_168839_168842_168885_168870_168848_168907_168851_476662_487594_19434_19415_19427_575809_366830_366838_620798_624001_498277_358596_358620_718737_667673_617828_733369_521067_321531_321537_321523_526723_68910_578655_182133_612303_174787_174813_174769_174807_174837_174819_562254_489394_69419_485733_211929_211891_211907_211960_211879_211918_211946_211935_211949_670145_600108_69475_69445_69457_69465_442326_89415_670254_26410_26437_515116_224754_224781_224732_224712_224758_224772_224684_554283_653488_223739_223735_223742_702006_735358_384310_384319_466060_176148_176106_176113_660323_19215_19200_19206_19223_433190_221741_221748_221716_221705_543628_650746_450382_224595_542366_680022_321612_677330_225129_225122_225104_225138_225117_225158_730499_722049_323887_323915_460209_312176_312185_312193_413621_319566_319557_319618_319598_319631_319604_319571_319607_658340_721511_361922_361946_361936_361912_532895_278911_278928_278944_490270_222337_222344_222358_222354_723792_368258_368265_521357_265628_694208_118962_118975_118932_118944_475305_10476_10469_10464_480489_369877_369882_394954_225090_225096_225106_225069_225042_225019_225081_225057_735858_311756_311793_311804_311776_311808_311801_539497_167849_725886_240888_240873_526618_26036_700959_396016_172638_172666_172628_172610_172619_172678_172657_603425_119241_119259_119285_119295_652308_69028_69015_451832_317983_317950_317945_317936_317975_317967_500406_274381_406735_20684_667510_120633_120621_120651_120636_120641_691247_747141_218422_218446_218430_218455_218413_636183_223795_223791_223819_223829_223825_223798_638501_331007_729226_222569_734190_170774_170765_170769_170728_170739_170747_656324_697635_759840_602140_87632_87647_87651_87667_87673_714108_594169_540510_562100_11278_11248_11283_11285_11267_753036_729640_378167_378206_378179_378216_378196_699288_372099_372076_372062_372088_372058_747491_219531_219523_219490_219497_219511_219537_617109_326681_326666_474876_143375_609592_573435_173905_173915_173909_531517_681874_229310_229318_229329_229301_514820_66453_66467_744186_71199_71205_71181_536100_371949_371932_539736_517542_65980_66001_66023_66010_66026_66015_693764_322103_521430_322129_547062_228324_228311_746002_178220_178153_178181_178166_178233_178194_178212_766707_219037_219009_218998_394700_184823_184827_184832_494344_127172_127104_127179_127130_127150_127146_127112_439135_477423_26724_608060_67252_67267_67240_67207_67260_67250_444850_719045_231305_231264_231296_231287_231276_632929_274809_274792_274806_531228_732918_328728_328735_635526_331670_331681_475631_63141_63153_63133_63131_630521_234510_702562_123048_122993_123031_123019_123004_729005_72748_473326_365724_365752_365702_365747_365706_365711_434468_144620_676892_335515_335501_485655_399545_219578_219584_555497_369655_369646_369651_369666_561797_234580_234592_234632_748523_175337_175354_175370_175382_175343_175364_175375_591145_442310_377209_377242_377233_377195_377240_377237_620721_584868_66607_767796_72716_509326_232099_232092_232071_232111_232105_232085_420102_413248_143943_143940_433831_174794_174792_174773_174800_174789_610442_21308_21303_748710_369887_411417_376954_376941_713378_296931_663427_457779_221192_221207_221199_221171_543826_321555_537013_483927_597759_77347_77378_77361_77340_77386_742467_178392_401112_171721_171673_171663_171691_171700_171728_171680_171697_758337_189241_507165_619031_183027_183035_183023_578618_593261_458157_22014_21976_22017_22009_22058_22047_22021_22039_22036_22031_639762_183764_183738_600653_64827_64955_64883_64837_64859_64896_64913_64871_64900_443127_328772_735029_275259_634493_220234_220285_220221_220210_220227_569003_62677_62691_62667_62660_62664_62698_408107_590162_74971_596905_271439_271450_271429_271394_271397_271455_271463_271434_271443_674641_69683_69673_69692_69664_69659_720681_420841_19237_19216_19240_19233_544268_400382_169059_169034_711168_370571_370589_537062_325180_650842_217268_217234_217254_217249_217237_217263_631440_227018_227054_227049_515245_279961_279975_279949_279967_602165_126066_126043_126050_126063_676120_120715_120744_479062_596243_271114_559039_170435_170440_170490_170501_170532_170470_621282_124572_540711_667645_68550_68543_68518_68528_68556_616159_15789_15772_15792_15781_15769_652034_12833_12821_12797_752160_229557_229579_229547_418930_220215_220233_220226_658209_182894_182888_182913_507263_325540_325554_325526_325523_325534_485611_738811_224034_222331_222316_222308_628649_66939_66927_66942_705602_179396_179402_179389_179424_466059_322596_322594_322584_716872_319532_492105_496285_732287_219976_219948_219997_219989_219963_219953_715699_719284_182170_182140_182151_182195_182176_668961_367874_367878_367887_367864_632097_747156_169642_169652_169665_169635_169630_474404_64560_64502_64524_64477_64535_64571_64565_64552_756381_436490_228589_228582_228601_228585_663766_223182_520997_164379_657822_211040_211046_211061_211056_520690_127625_510156_760064_16807_16856_16811_16844_499115_318734_318748_318744_318730_318755_318689_318797_318788_318813_318770_655332_505749_12990_13003_12986_12972_690739_126751_126756_663860_116494_116505_116478_116529_116498_116488_704564_15267_698233_324684_589902_468720_280680_667765_273089_273096_273062_273053_273072_692837_275709_519714_375554_593939_654810_133088_133095_133075_418709_59715_59734_59722_59781_59778_59739_59768_710539_224104_399544_661046_313081_313101_313095_313084_484731_364078_364093_364040_364064_364057_364071_603311_36661_36651_612334_567213_220697_220683_220713_644189_76669_610916_68452_68443_68410_68432_524262_22801_22786_651379_396438_23261_23247_408013_282654_282636_282648_282671_706318_131923_131919_553081_81293_511091_288704_666866_24192_24190_24230_24215_24185_667291_29235_29238_29215_29229_767653_556856_174950_750631_648486_716411_762688_455120_420404_231456_231475_231480_231463_231405_231440_231447_231471_231409_490056_618031_594141_285203_285223_285181_285186_285196_285209_285215_401211_333799_333813_333805_333782_750459_488828_716860_369418_369443_369423_369400_369411_369437_403376_401960_661269_23824_23838_23791_23817_23833_23829_695408_29369_29391_29394_29383_541419_665551_168824_168818_168828_168843_168847_168838_168831_168813_709894_138835_627008_88371_701462_288359_516969_174702_698603_325840_325857_325836_325826_325875_721970_171309_171289_171292_626362_325381_325368_325370_325385_325386_325377_452795_179764_179786_179753_179726_702558_66733_66710_66722_66718_66714_66735_66727_521290_288231_288225_540315_572442_19721_19705_546260_370301_370317_370306_744894_216957_216961_397396_412505_218020_218042_218016_218027_437890_120780_674235_324089_324086_577606_276574_276596_276592_676344_664426_173549_582455_508724_78104_78121_78107_78101_78082_78072_402407_134076_134095_545711_333041_333023_563866_13943_13960_13952_13940_13970_13947_419923_368494_368502_368489_368471_368447_487753_446733_22254_22241_22265_634130_221180_699194_118411_118414_118432_416673_418632_265019_265011_265023_265001_265004_472896_331750_478863_67119_67130_671480_627590_63329_63367_63375_63386_63335_409552_169937_169953_169968_169976_169896_169902_497782_405201_24452_24474_24487_511813_399332_268141_268176_268152_268163_268145_469365_661638_8842_8870_8861_8879_542704_162911_162889_162880_162883_697703_264791_264813_264806_561635_274136_274132_274180_274173_550526_459260_369868_369860_369888_369846_369865_622762_271230_271243_401481_681584_464269_371899_371903_371910_685344_330499_674478_178270_503415_225017_225009_225004_579409_243795_642846_280562_280571_580358_121742_742686_322764_322779_488810_384991_385004_551279_187625_187611_187620_187593_458928_444803_399241_82386_82332_82339_82342_82363_698520_127784_127743_127754_127789_127734_502475_376591_376548_376577_376573_376562_529737_638185_280450_426706_126687_126727_126714_126670_126719_435727_326032_326026_326053_326016_468671_589534_374570_374585_671354_374993_374996_455294_169514_169498_169542_169529_169546_540884_363280_363255_363273_670026_272302_272318_272328_697756_705354_483004_330673_330658_330653_330714_330756_330697_330749_541133_536478_122414_122408_122397_402653_685249_581713_343297_343305_600941_91904_403002_224258_224251_224267_224247_224235_760155_226396_226414_412250_218385_218373_410611_169776_169791_169769_169761_169750_169786_576708_596788_57332_57291_57323_57318_619153_612097_239834_509012_705376_273202_470488_66833_66796_66809_66818_66842_66825_751990_173368_173335_173361_517388_218876_218883_218868_218894_218906_218903_701073_169323_169312_533362_173115_173098_173119_173106_536578_673497_129035_129051_129086_129080_129075_129020_129028_129043_471227_596915_651942_279631_279640_279668_279690_688326_286028_286013_517472_382677_402595_76980_541134_180215_180272_180250_180230_705086_514892_384001_384020_646551_377078_377127_377134_377137_377097_377089_377120_377084_519972_284482_759382_137351_137361_137372_137354_137345_732224_371083_371091_371060_371114_371056_371122_702605_518870_141058_141049_141052_685962_346040_347923_347938_347932_346047_347911_594256_133181_133243_133219_133200_133230_457584_476046_234114_234053_234073_234088_234099_234043_234117_234110_580877_328572_328567_607539_625296_285493_285461_724868_278894_278881_278887_576458_137664_137652_457192_239521_239530_239538_475716_731867_644597_131543_131513_131510_131518_131526_597325_28002_461483_279692_279699_427607_400843_220124_220127_220113_220118_220104_764859_29093_29096_29112_29076_29048_29130_29054_557212_515858_15605_15582_15654_15619_15637_15626_15647_457072_180565_560284_282044_282058_410057_671216_317524_317545_317515_317518_317522_608376_15888_15873_15861_15838_15883_15897_15856_519152_125646_763112_18766_754170_119653_119659_764717_125505_125514_125452_125479_125457_125522_531338_220981_220972_220973_220958_220909_220947_220931_220911_220966_763904_684539_367040_661505_323865_323852_323846_702419_218614_218604_218625_218612_218636_668545_720792_737071_69258_69246_69267_646090_234082_234125_234143_234101_234086_234076_234152_234119_646911_279120_279124_427826_67020_67026_67017_558520_557539_717337_273865_273843_273833_273824_554779_677743_172813_172806_438728_276837_726486_173183_173210_173188_173230_173180_708554_327776_327814_327792_327762_327799_432556_324825_764043_70524_70507_70543_70498_70528_70548_70513_663719_632443_272123_272119_272038_272052_272103_272059_272096_491700_12654_12676_12665_12640_12648_12658_12632_607849_526424_369040_369052_369077_369046_692525_178818_746597_680849_271825_271776_271800_271790_271793_703938_17937_17948_17909_17932_17943_17881_17929_17897_17888_17911_500233_20148_20126_20171_20141_660310_327743_327728_574759_372994_373030_373035_373008_373014_373041_479885_324853_324871_324899_324888_324864_502362_516066_81113_81117_420165_277357_277371_277331_277342_478558_72614_72619_526418_524382_180879_744468_22253_22221_22250_22238_568186_16183_16155_16145_16129_16172_16138_16187_656053_20421_20418_20438_20411_20404_20415_672365_74317_74285_74282_74268_74308_760112_372747_372761_372716_372742_372753_372721_372729_482063_291119_291133_291115_622036_509759_136099_136175_136125_136103_136119_136190_136110_136149_136185_697584_447574_176962_176944_176954_402453_127160_127166_127162_127184_695177_75590_75600_520022_74211_74193_652800_278991_278985_278958_278971_279000_505130_580249_134784_134772_134768_134756_134799_660997_229727_229739_229733_229712_464220_509734_616269_64876_64862_64872_601787_124453_124465_124446_124470_419988_172251_172270_172245_172238_438729_323635_323615_496522_14967_14952_671588_120023_120008_120057_120051_120031_603963_180210_180208_536108_678850_68290_68280_439798_331111_331119_331107_331105_486442_754678_162501_162469_162472_162482_162478_162467_715620_59374_59386_59340_59403_59345_59365_673517_624774_232002_623580_368394_368381_368374_368390_599974_14251_14291_14268_14324_14314_14306_14256_400455_174973_174981_566961_18134_18127_702479_69758_69766_69737_69754_69748_69740_69806_575587_439797_275993_275998_275954_275966_428946_598942_127956_127961_127937_127949_127977_127926_397055_492954_174665_500929_218862_218855_575519_30320_30299_30332_30349_558126_140379_140366_407714_67898_67890_457947_317943_317977_317962_317968_461689_367570_367537_496979_125953_125938_125961_125948_555093_221089_221124_221211_221141_221176_221101_221162_221189_221115_221133_221067_754118_7902_7908_421882_313335_313325_313321_674543_10129_10108_10117_10102_558162_57864_57878_57866_57913_57850_57869_57897_606012_460812_282412_282396_594466_671461_384666_743531_314194_314216_314201_314176_314169_489942_324224_324255_324244_324229_324237_446659_226629_226646_226672_226666_226662_226637_226689_226632_226681_600803_275074_275067_275039_275051_692162_278208_501301_236890_236898_236904_236887_735846_453293_319535_463966_127350_127337_732946_319869_319875_319855_319836_319852_319882_319862_319830_434586_486495_593177_31534_31502_31525_31514_31543_723399_236416_482412_333237_464524_372857_372862_372868_626683_328172_328186_328196_328153_328144_328178_328141_725740_609915_700068_174106_174176_174171_174162_174136_761632_372312_372319_699114_224418_224462_224437_224467_224453_224411_224394_488267_761340_532050_405257_457449_273681_273666_273696_470036_275206_275221_275187_275200_744227_67509_67424_67429_67500_67458_67489_67477_635748_28512_571739_20860_20852_20865_20856_613289_188280_722645_124872_124900_124896_124860_453363_416448_740494_375357_375336_375366_375360_375379_590847_218711_218684_218702_218729_218706_615614_36364_580380_283958_484726_177247_177250_399996_181352_181384_543600_325193_325187_325177_325202_325171_736607_225008_225014_545087_66569_66586_66603_66592_677284_80053_80066_650718_223624_223650_223629_223607_223634_547167_446563_175971_175994_175997_176023_176034_176062_176040_543775_120914_120899_120842_120878_120871_430376_29338_29341_560804_365590_572993_441852_715166_124111_124169_124123_124159_124141_124108_497922_16771_584746_331586_331593_331601_549325_181832_181812_181838_689822_11692_11723_11707_11738_11743_11716_11704_746849_65236_761577_70759_70765_526729_540277_25422_25462_25447_25471_25434_25441_445000_132412_132406_132451_132423_132431_132440_132461_620449_598779_130488_130480_431491_220189_220197_220164_220170_500337_120347_120334_646355_374671_374638_374646_374663_374655_397486_30650_30667_30684_30678_429214_284660_706869_134617_134626_134656_134633_134673_134644_726489_369625_369636_369630_516715_324998_325024_325064_325015_325027_652428_322604_322609_322619_322613_322606_482455_606807_321969_321959_321948_321979_541610_322311_322329_322289_322322_322273_477623_370145_715685_217906_217882_217949_217912_217895_217933_662423_448073_14840_14816_14844_14768_14852_14777_14876_14834_14772_14804_14879_728765_322254_322247_322302_322296_322271_322278_642952_70136_70148_70159_615513_497913_67921_67927_67953_67935_67919_588094_690257_512052_68707_68669_68704_68679_68699_631053_174017_174061_174022_174044_174004_174047_174056_174027_758269_62879_62926_62902_62875_62887_658124_328433_328429_328426_328419_328423_517782_328330_720337_227624_227619_447184_223281_223319_223306_223293_223310_592491_65573_439342_180209_180187_664521_331570_331565_540720_76460_76468_76438_743923_370271_528386_76146_76155_76161_76172_468965_74310_561457_271529_659355_131229_131261_131248_131257_501275_330752_330778_330761_330758_487434_412349_279862_609681_24002_23992_24005_501564_334523_334515_463282_273124_273116_647612_177954_539621_20228_20222_20210_20203_20197_20192_572540_227749_227752_227757_537621_416469_364435_364473_364414_364433_364479_364456_364429_364463_364484_364440_653652_658000_551215_165836_165824_165820_165840_507824_215499_215494_215487_215501_760777_7712_6445_7760_6455_7773_6449_7765_7757_7735_588109_271055_271048_596354_121943_578006_162857_162843_162862_162866_697032_712129_176262_176289_176274_474895_721997_673676_25785_25781_25767_25762_497088_436172_80565_80558_756387_63003_62996_62988_63006_62970_652240_222574_222550_418060_219968_220001_219977_220028_219982_220024_764659_181029_755784_78159_78169_78198_78162_495814_13668_13678_13661_13664_13683_621420_228950_573872_180785_180794_180777_180788_180803_180759_509929_136480_136473_136470_437993_16231_16207_16197_671415_370751_370730_667505_698848_18061_18075_18079_18042_18046_716854_519079_78654_78673_78621_78662_78612_78638_78592_78599_530104_377406_437140_234077_234070_542512_369050_369042_369036_369059_369014_369072_448360_321277_710373_627679_627617_719478_644617_129105_129054_129081_129098_129076_591820_274742_274770_274757_274753_274732_274728_481588_331353_331388_331362_578650_132465_625795_329018_329032_329009_329023_329014_329053_442865_327574_327545_327567_617015_668345_220194_220174_220148_220156_576020_629500_139572_139580_469441_29182_29104_29108_29138_29094_29111_29132_465162_332249_332255_586443_279277_620654_332822_332787_332814_553766_319924_319897_319880_319853_319908_319890_319864_594901_89447_89460_89427_565617_25929_25924_25940_25898_25890_25957_25907_25885_25967_549406_35996_35978_35970_35986_702103_74741_74695_74724_74684_74710_475811_21458_440675")
        // {
        //     continue;
        // }
        if (length_of_contigs[contig_name] > length_of_longest_read){

            //list all paths from the left end of the contig
            vector<vector<pair<string,char>>> paths_left = list_all_paths_from_this_end(contig_name, 0, (int) length_of_longest_read/2, linked, length_of_contigs);
            
            //check if all paths lead to the same contig
            pair<string,char> contig_to_reach = {"", 0};
            bool all_paths_to_same_contig = false;
            for (auto p: paths_left){
                if (p[0].first != contig_name && contig_to_reach.first == ""){
                    contig_to_reach = p[0];
                    all_paths_to_same_contig = true;
                }
                else if (p[0].first != contig_name && contig_to_reach != p[0]){
                    all_paths_to_same_contig = false;
                }
            }

            //if there is all paths lead there check to see if it is reciprocal
            if (all_paths_to_same_contig){
                vector<vector<pair<string,char>>> paths_right = list_all_paths_from_this_end(contig_to_reach.first, 1-contig_to_reach.second, (int) length_of_longest_read/2, linked, length_of_contigs);
                bool all_paths_reciprocal = true;
                for (auto p: paths_right){
                    if (p[0].first != contig_to_reach.first &&  (p[0].first != contig_name || p[0].second != 1)){
                        all_paths_reciprocal = false;
                    }
                }

                if (all_paths_reciprocal){

                    // cout << "popping bubble between " << contig_name << " and " << contig_to_reach.first << "\n";
                    
                    //now compute the average coverage of each path of the bubble
                    vector<float> coverages (paths_left.size(), 0);
                    int best_path = 0;
                    float best_coverage = 0;
                    for (int i = 0 ; i < paths_left.size() ; i++){
                        int length_of_path = 0;
                        for (auto p: paths_left[i]){
                            coverages[i] += coverage[p.first]*length_of_contigs[p.first];
                            length_of_path += length_of_contigs[p.first];
                        }
                        coverages[i] /= length_of_path;
                        if (coverages[i] >= best_coverage && paths_left[i][0].first != contig_name){
                            best_coverage = coverages[i];
                            best_path = i;
                        }
                    }

                    //take the path with the highest coverage, go through it and detach all the other paths
                    std::set<pair<pair<string,char>,pair<string,char>>> links_to_keep;
                    for (int co = 1 ; co < paths_left[best_path].size() ; co++){
                        links_to_keep.insert({{paths_left[best_path][co-1].first, 1-paths_left[best_path][co-1].second}, paths_left[best_path][co]});
                    }
                    for (int i = 0 ; i < paths_left.size() ; i++){
                        if (i != best_path){
                            for (int co = 1 ; co < paths_left[i].size() ; co++){
                                //check if it is not in the links to keep
                                if (links_to_keep.find({{paths_left[i][co-1].first, 1-paths_left[i][co-1].second}, paths_left[i][co]}) == links_to_keep.end() 
                                    && links_to_keep.find({paths_left[i][co],{paths_left[i][co-1].first, 1-paths_left[i][co-1].second} }) == links_to_keep.end()  ){
                                    links_to_delete.insert({paths_left[i][co],{paths_left[i][co-1].first, 1-paths_left[i][co-1].second} });
                                    links_to_delete.insert({{paths_left[i][co-1].first, 1-paths_left[i][co-1].second}, paths_left[i][co]});
                                }
                            }
                        }
                    }
                }
            }

            //do the same thing for the right end of the contig
            paths_left = list_all_paths_from_this_end(contig_name, 1, (int) length_of_longest_read/2, linked, length_of_contigs);
            all_paths_to_same_contig = false;
            contig_to_reach = {"", 0};
            for (auto p: paths_left){
                if (contig_to_reach.first == ""){
                    contig_to_reach = p[0];
                    all_paths_to_same_contig = true;
                }
                else if (contig_to_reach != p[0]){
                    all_paths_to_same_contig = false;
                }
            }

            // cout << "here are all the paths from the right end of " << contig_name.substr(0,30) << "\n";
            // for (auto p: paths_left){
            //     for (auto c: p){
            //         cout << c.first.substr(0,30)  << " ";
            //     }
            //     cout << "\n";
            // }

            //if there is all paths lead there check to see if it is reciprocal
            if (all_paths_to_same_contig){
                vector<vector<pair<string,char>>> paths_right = list_all_paths_from_this_end(contig_to_reach.first, 1-contig_to_reach.second, (int) length_of_longest_read/2, linked, length_of_contigs);
                bool all_paths_reciprocal = true;
                for (auto p: paths_right){
                    if (p[0].first != contig_name || p[0].second != 0){
                        all_paths_reciprocal = false;
                    }
                }

                if (all_paths_reciprocal){

                    // cout << "popping bubble between " << contig_name << " and " << contig_to_reach.first << "\n";

                    //now compute the average coverage of each path of the bubble
                    vector<float> coverages (paths_left.size(), 0);
                    int best_path = 0;
                    float best_coverage = 0;
                    for (int i = 0 ; i < paths_left.size() ; i++){
                        int length_of_path = 0;
                        for (auto p: paths_left[i]){
                            coverages[i] += coverage[p.first]*length_of_contigs[p.first];
                            length_of_path += length_of_contigs[p.first];
                        }
                        coverages[i] /= length_of_path;
                        if (coverages[i] >= best_coverage){
                            best_coverage = coverages[i];
                            best_path = i;
                        }
                    }

                    //take the path with the highest coverage, go through it and detach all the other paths
                    std::set<pair<pair<string,char>,pair<string,char>>> links_to_keep;
                    for (int co = 1 ; co < paths_left[best_path].size() ; co++){
                        links_to_keep.insert({{paths_left[best_path][co-1].first, 1-paths_left[best_path][co-1].second}, paths_left[best_path][co]});
                    }
                    for (int i = 0 ; i < paths_left.size() ; i++){
                        if (i != best_path){
                            for (int co = 1 ; co < paths_left[i].size() ; co++){
                                //check if it is not in the links to keep
                                if (links_to_keep.find({{paths_left[i][co-1].first, 1-paths_left[i][co-1].second}, paths_left[i][co]}) == links_to_keep.end() 
                                    && links_to_keep.find({paths_left[i][co],{paths_left[i][co-1].first, 1-paths_left[i][co-1].second} }) == links_to_keep.end()  ){
                                    links_to_delete.insert({paths_left[i][co],{paths_left[i][co-1].first, 1-paths_left[i][co-1].second} });
                                    links_to_delete.insert({{paths_left[i][co-1].first, 1-paths_left[i][co-1].second}, paths_left[i][co]});
                                }
                            }
                        }
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

void load_GFA(string gfa_file, vector<Segment> &segments, unordered_map<string, int> &segment_IDs){
    //load the segments from the GFA file
    
    //in a first pass index all the segments by their name
    ifstream gfa(gfa_file);
    string line;
    while (getline(gfa, line)){
        if (line[0] == 'S'){
            long int pos_in_file = (long int) gfa.tellg() - line.size() - 1;
            stringstream ss(line);
            string nothing, name;
            ss >> nothing >> name;

            double coverage = 0;
            //try to find a DP: tag
            string tag;
            while (ss >> tag){
                if (tag.substr(0, 3) == "DP:" || tag.substr(0, 3) == "km:"){
                    coverage = std::stof(tag.substr(5, tag.size() - 5));
                }
            }

            Segment s(name, segments.size(), vector<pair<vector<pair<int,int>>, vector<string>>>(2), pos_in_file, coverage);

            segment_IDs[name] = s.ID;
            segments.push_back(s);
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

void merge_adjacent_contigs(vector<Segment> &old_segments, vector<Segment> &new_segments, string original_gfa_file, bool rename){

    set<int> already_looked_at_segments; //old IDs of segments that have already been looked at and merged (don't want to merge them twice)
    int seg_idx = 0;
    unordered_map<pair<int,int>,pair<int,int>> old_ID_to_new_ID; //associates (old_id, old end) with (new_id_new_end)
    int number_of_merged_contigs = 0;
    set<pair<pair<pair<int,int>, pair<int,int>>,string>> links_to_add; //list of links to add, all in old IDs
    for (Segment old_seg : old_segments){

        // if (old_seg.name != "1665420_0"){
        //     seg_idx++;
        //     continue;
        // }

        if (already_looked_at_segments.find(old_seg.ID) != already_looked_at_segments.end()){
            seg_idx++;
            continue;
        }
        //check if it has either at least two neighbors left or that its neighbor left has at least two neighbors right
        cout << "in merge, looking at segment " << seg_idx << " out of " << old_segments.size() << "\r" << std::flush;
        bool dead_end_left = false;
        if (old_seg.links[0].first.size() != 1 || old_segments[old_seg.links[0].first[0].first].links[old_seg.links[0].first[0].second].first.size() != 1 || old_segments[old_seg.links[0].first[0].first].ID == old_seg.ID){
            dead_end_left = true;
        }

        bool dead_end_right = false;
        if (old_seg.links[1].first.size() != 1 || old_segments[old_seg.links[1].first[0].first].links[old_seg.links[1].first[0].second].first.size() != 1 || old_segments[old_seg.links[1].first[0].first].ID == old_seg.ID){
            dead_end_right = true;
        }

        if (!dead_end_left && !dead_end_right){ //means this contig is in the middle of a long haploid contig, no need to merge
            seg_idx++;
            continue;
        }


        //create a new contig
        if (dead_end_left && dead_end_right){
            already_looked_at_segments.insert(old_seg.ID);
            string name = old_seg.name;
            if (rename){
                name = std::to_string(number_of_merged_contigs);
                number_of_merged_contigs++;
            }
            new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), old_seg.get_coverage()));
            old_ID_to_new_ID[{old_seg.ID, 0}] = {new_segments.size() - 1, 0};
            old_ID_to_new_ID[{old_seg.ID, 1}] = {new_segments.size() - 1, 1};
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
            // cout << "yupee" << endl;
        }
        else if (dead_end_left){
            
            //let's see how far we can go right
            vector<string> all_names = {old_seg.name};
            vector<string> all_seqs = {old_seg.get_seq(original_gfa_file)};
            vector<double> all_coverages = {old_seg.get_coverage()};
            int current_ID = old_seg.ID;
            int current_end = 1;

            while (old_segments[current_ID].links[current_end].first.size() == 1 && old_segments[old_segments[current_ID].links[current_end].first[0].first].links[old_segments[current_ID].links[current_end].first[0].second].first.size() == 1){
                already_looked_at_segments.insert(current_ID);
                string cigar = old_segments[current_ID].links[current_end].second[0];
                int tmp_current_end = 1-old_segments[current_ID].links[current_end].first[0].second;
                current_ID = old_segments[current_ID].links[current_end].first[0].first;
                current_end = tmp_current_end;
                all_names.push_back(old_segments[current_ID].name);
                string seq = old_segments[current_ID].get_seq(original_gfa_file);
                //now reverse complement if current_end is 1
                if (current_end == 0){
                    seq = reverse_complement(seq);
                }
                //trim the sequence if there is a CIGAR
                int num_matches = std::stoi(cigar.substr(0, cigar.find_first_of("M")));
                all_seqs.push_back(seq.substr(num_matches, seq.size()-num_matches));
                all_coverages.push_back(old_segments[current_ID].get_coverage());
            }
            already_looked_at_segments.insert(current_ID);

            //create the new contig
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
            for (double coverage : all_coverages){
                new_coverage += coverage;
            }
            new_coverage = new_coverage/all_coverages.size();
            string name = new_name;
            if (rename){
                name = std::to_string(number_of_merged_contigs);
                number_of_merged_contigs++;
            }
            new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), new_coverage));
            new_segments[new_segments.size()-1].seq = new_seq;

            //add the links
            int idx_link = 0;
            for (pair<int,int> link : old_seg.links[0].first){
                links_to_add.insert({{{old_seg.ID, 0}, link}, old_seg.links[0].second[idx_link]});
                idx_link++;
            }
            idx_link = 0;
            for (pair<int,int> link : old_segments[current_ID].links[current_end].first){
                links_to_add.insert({{{current_ID, current_end}, link}, old_segments[current_ID].links[current_end].second[idx_link]});
                idx_link++;
            }

            old_ID_to_new_ID[{old_seg.ID, 0}] = {new_segments.size() - 1, 0};
            old_ID_to_new_ID[{current_ID, current_end}] = {new_segments.size() - 1, 1};

            // cout << "hurray" << endl;
        }
        else{
            
            //let's see how far we can go left
            vector<string> all_names = {old_seg.name};
            string seq = old_seg.get_seq(original_gfa_file);
            vector<string> all_seqs = {reverse_complement(seq)};
            vector<double> all_coverages = {old_seg.get_coverage()};
            int current_ID = old_seg.ID;
            int current_end = 0;

            // cout << "exploring all the contigs left" << endl;
            // cout << "first exploring the link between " << old_segments[current_ID].name << " and " << old_segments[old_segments[current_ID].links[current_end].first[0].first].name << endl;
            
            while (old_segments[current_ID].links[current_end].first.size() == 1 && old_segments[old_segments[current_ID].links[current_end].first[0].first].links[old_segments[current_ID].links[current_end].first[0].second].first.size() == 1){
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

                // cout << "exploring the link between " << old_segments[current_ID].name << " and " << old_segments[old_segments[current_ID].links[current_end].first[0].first].name  << " " << old_segments[current_ID].links[current_end].first.size() << endl;
                // cout << "here are all the links : ";
                // for (pair<int,int> link : old_segments[current_ID].links[current_end].first){
                //     cout << link.first << " " << old_segments[link.first].name << " " << link.second << " ; ";
                // }
                // cout << endl;
            }
            already_looked_at_segments.insert(current_ID);

            //create the new contig
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
            for (double coverage : all_coverages){
                new_coverage += coverage;
            }
            new_coverage = new_coverage/all_coverages.size();
            string name = new_name;
            if (rename){
                name = std::to_string(number_of_merged_contigs);
                number_of_merged_contigs++;
            }
            new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), new_coverage));
            new_segments[new_segments.size()-1].seq = new_seq;

            //add the links
            int idx_link = 0;
            for (pair<int,int> link : old_seg.links[1].first){
                links_to_add.insert({{{old_seg.ID, 1}, link}, old_seg.links[1].second[idx_link]});
                idx_link++;
            }
            idx_link = 0;
            for (pair<int,int> link : old_segments[current_ID].links[current_end].first){
                links_to_add.insert({{{current_ID, current_end}, link}, old_segments[current_ID].links[current_end].second[idx_link]});
                idx_link++;
            }

            old_ID_to_new_ID[{old_seg.ID, 1}] = {new_segments.size() - 1, 0};
            old_ID_to_new_ID[{current_ID, current_end}] = {new_segments.size() - 1, 1};

            // cout << "yay" << endl;
        }

        seg_idx++;
    }

    //some contigs are left: the ones that were in circular rings... go through them and add them
    for (Segment old_seg : old_segments){
        if (already_looked_at_segments.find(old_seg.ID) == already_looked_at_segments.end()){
            int current_ID = old_seg.ID;
            int current_end = 1;
            vector<string> all_names = {old_seg.name};
            vector<string> all_seqs = {old_seg.get_seq(original_gfa_file)};
            vector<double> all_coverages = {old_seg.get_coverage()};
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
                for (double coverage : all_coverages){
                    new_coverage += coverage;
                }
                new_coverage = new_coverage/all_coverages.size();
                string name = new_name;
                if (rename){
                    name = std::to_string(number_of_merged_contigs);
                    number_of_merged_contigs++;
                }
                new_segments.push_back(Segment(name, new_segments.size(), old_seg.get_pos_in_file(), new_coverage));
                new_segments[new_segments.size()-1].seq = new_seq;

                //add the links
                int idx_link = 0;
                for (pair<int,int> link : old_seg.links[0].first){
                    links_to_add.insert({{{old_seg.ID, 0}, link}, old_seg.links[0].second[idx_link]});
                    idx_link++;
                }

                old_ID_to_new_ID[{old_seg.ID, 1}] = {new_segments.size() - 1, 0};
                old_ID_to_new_ID[{current_ID, current_end}] = {new_segments.size() - 1, 1};
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
        gfa << "S\t" << s.name << "\t" << s.get_seq(gfa_input) << "\tDP:f:" << s.get_coverage() <<  "\n";
    }
    for (Segment s : segments){
        for (int end = 0 ; end < 2 ; end++){
            for (int neigh = 0 ; neigh < s.links[end].first.size() ; neigh++){

                //to make sure the link is not outputted twice
                if (s.ID > s.links[end].first[neigh].first || (s.ID == s.links[end].first[neigh].first && end > s.links[end].first[neigh].second) ){
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


