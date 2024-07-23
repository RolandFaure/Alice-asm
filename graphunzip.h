#ifndef GRAPHUNZIP_H
#define GRAPHUNZIP_H

#include <iostream>
#include <vector>
#include <string>
#include <mutex>
#include <set>
#include <string>
#include <fstream>
#include <sstream>

#include "robin_hood.h"
#include "basic_graph_manipulation.h"

// namespace std {
//     template <>
//     class hash<std::pair<int, int>> {
//     public:
//         size_t operator()(const std::pair<int, bool>& pair) const {
//             return std::hash<int>()(pair.first) ^ std::hash<int>()(pair.second);
//         }
//     };
// }


// class Segment{
//     public:

//         std::string name;
//         int ID;
//         std::vector<std::pair<std::vector<std::pair<int,int>>, std::vector<std::string>>>  links; //first vector is for the links to the left, second vector is for the links to the right. Each link is the ID of the neighbor and its end (0 for left, 1 for right) and the CIGAR

//         Segment(){};

//         Segment(std::string name, int ID, long int pos_in_file, double coverage){
//             this->name = name;
//             this->ID = ID;
//             this->pos_in_file = pos_in_file;
//             this->coverage = coverage;
//             this->haploid = false;
//             this->links = std::vector<std::pair<std::vector<std::pair<int,int>>, std::vector<std::string>>>(2);
//         }

//         Segment(std::string name, int ID, std::vector<std::pair<std::vector<std::pair<int,int>>, std::vector<std::string>>> links, long int pos_in_file, double coverage){
//             this->name = name;
//             this->ID = ID;
//             this->links = links;
//             this->pos_in_file = pos_in_file;
//             this->coverage = coverage;
//             this->haploid = false;
//             this->links = std::vector<std::pair<std::vector<std::pair<int,int>>, std::vector<std::string>>>(2);
//         }

//         bool is_haploid(){return this->haploid;}
//         std::vector<std::pair<int,bool>> get_consensus_left(){return this->consensus_left;}
//         std::vector<std::pair<int,bool>> get_consensus_right(){return this->consensus_right;}
//         long int get_pos_in_file(){return this->pos_in_file;}
//         double get_coverage(){return this->coverage;}

//         std::string get_seq(std::string& gfa_file){
            
//             if (seq != ""){
//                 // cout << "the seq is already loaded" << endl;
//                 return seq;
//             }

//             std::ifstream gfa(gfa_file);
//             gfa.seekg(pos_in_file);
//             std::string line;
//             std::getline(gfa, line);
//             std::stringstream ss(line);
//             std::string nothing, name, seq;
//             ss >> nothing >> name >> seq;
//             gfa.close();
//             return seq;
//         }

//         void add_neighbor(std::vector<std::pair<int,bool>> new_neighbor, bool left);

//         void decrease_coverage(double coverage_out){
//             coverage -= coverage_out;
//             if (coverage < 1){
//                 coverage = 1;
//             }
//         }

//         void compute_consensuses();

//         std::vector<std::vector<std::pair<int,bool>>> get_strong_neighbors_left();
//         std::vector<std::vector<std::pair<int,bool>>> get_strong_neighbors_right();

//         std::vector<std::pair<std::vector<std::pair<int,bool>>, int>> get_neighbors_left(){return neighbors_left;}
//         std::vector<std::pair<std::vector<std::pair<int,bool>>, int>> get_neighbors_right(){return neighbors_right;}

//         //the hash of the segment is the hash of the name
//         size_t hash() const{
//             return std::hash<std::string>{}(name);
//         }

//         std::string seq;

//     private:
//         //consensus sequences of contigs left and right
//         std::vector<std::pair<int,bool>> consensus_left;
//         std::vector<std::pair<int,bool>> consensus_right;

//         std::vector<std::vector<std::pair<int,bool>>> representative_neighbors_left;
//         std::vector<std::vector<std::pair<int,bool>>> representative_neighbors_right;

//         std::vector<std::pair<std::vector<std::pair<int,bool>>, int>> neighbors_left; //each path is a std::vector of std::pairs (ID, orientation) of the contigs in the path and a strength (number of reads)
//         std::vector<std::pair<std::vector<std::pair<int,bool>>, int>> neighbors_right; //each path is a std::vector of std::pairs (ID, orientation) of the contigs in the path and a strength (number of reads)

//         long int pos_in_file;
//         double coverage;

//         bool haploid;

// };

void load_GFA(std::string gfa_file, std::vector<Segment> &segments, robin_hood::unordered_map<std::string, int> &segment_IDs);
void merge_adjacent_contigs(std::vector<Segment> &old_segments, std::vector<Segment> &new_segments, std::string original_gfa_file, bool rename);
void output_graph(std::string gfa_output, std::string gfa_input, std::vector<Segment> &segments);


#endif