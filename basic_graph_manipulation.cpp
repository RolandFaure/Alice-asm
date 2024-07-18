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
        if (length_of_contigs[contig_name] > length_of_longest_read){

            //list all paths from the left end of the contig
            vector<vector<pair<string,char>>> paths_left = list_all_paths_from_this_end(contig_name, 0, (int) length_of_longest_read/2, linked, length_of_contigs);
            
            //check if all paths lead to the same contig
            pair<string,char> contig_to_reach = {"", 0};
            bool all_paths_to_same_contig = false;
            for (auto p: paths_left){
                if (contig_to_reach.first == ""){
                    contig_to_reach = p[0];
                    all_paths_to_same_contig = true;
                }
                else if (contig_to_reach != p[0]){
                    all_paths_to_same_contig = false;
                }
            }

            //if there is all paths lead there check to see if it is reciprocal
            if (all_paths_to_same_contig){
                vector<vector<pair<string,char>>> paths_right = list_all_paths_from_this_end(contig_to_reach.first, 1-contig_to_reach.second, (int) length_of_longest_read/2, linked, length_of_contigs);
                bool all_paths_reciprocal = true;
                for (auto p: paths_right){
                    if (p[0].first != contig_name || p[0].second != 1){
                        all_paths_reciprocal = false;
                    }
                }

                if (all_paths_reciprocal){

                    cout << "popping bubble between " << contig_name << " and " << contig_to_reach.first << "\n";
                    
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

            //do the same thing for the right end of the contig
            paths_left = list_all_paths_from_this_end(contig_name, 1, (int) length_of_longest_read/2, linked, length_of_contigs);
            all_paths_to_same_contig = false;
            for (auto p: paths_left){
                if (contig_to_reach.first == ""){
                    contig_to_reach = p[0];
                    all_paths_to_same_contig = true;
                }
                else if (contig_to_reach != p[0]){
                    all_paths_to_same_contig = false;
                }
            }

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

                    cout << "popping bubble between " << contig_name << " and " << contig_to_reach.first << "\n";

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
                if (tag.substr(0, 3) == "DP:"){
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

void merge_adjacent_contigs(vector<Segment> &old_segments, vector<Segment> &new_segments, string original_gfa_file){

    set<int> already_looked_at_segments; //old IDs of segments that have already been looked at and merged (don't want to merge them twice)
    int seg_idx = 0;
    unordered_map<pair<int,int>,pair<int,int>> old_ID_to_new_ID; //associates (old_id, old end) with (new_id_new_end)
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

        // cout << "dead end left: " << dead_end_left << " dead end right: " << dead_end_right << endl;

        if (!dead_end_left && !dead_end_right){ //means this contig is in the middle of a long haploid contig, no need to merge
            seg_idx++;
            continue;
        }

        //create a new contig
        if (dead_end_left && dead_end_right){
            already_looked_at_segments.insert(old_seg.ID);
            new_segments.push_back(Segment(old_seg.name, new_segments.size(), old_seg.get_pos_in_file(), old_seg.get_coverage()));
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
            new_segments.push_back(Segment(new_name, new_segments.size(), old_seg.get_pos_in_file(), new_coverage));
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
            new_segments.push_back(Segment(new_name, new_segments.size(), old_seg.get_pos_in_file(), new_coverage));
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


