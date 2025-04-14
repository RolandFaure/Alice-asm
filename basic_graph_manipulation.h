#ifndef BGM_H
#define BGM_H

#include <string>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include "robin_hood.h"
#include "rolling_hash.h"

std::string reverse_complement(std::string& seq);
void gfa_to_fasta(std::string gfa, std::string fasta);
void sort_GFA(std::string gfa);

void shave(std::string input_file, std::string output_file, int max_length);
void pop_and_shave_graph(std::string gfa_in, int abundance_min, int min_length, bool contiguity, int k, std::string gfa_out, int extra_coverage, int num_threads, bool single_genome);
void pop_bubbles(std::string gfa_in, int length_of_longest_read, std::string gfa_out);
void trim_tips_and_isolated_contigs(std::string gfa_in, int min_coverage, int min_length, std::string gfa_out);
void merge_adjacent_contigs_BCALM(std::string gfa_in, std::string gfa_out, int k, std::string path_to_bcalm, std::string path_convertToGFA, std::string path_tmp_folder);

void create_gaf_from_unitig_graph(std::string unitig_graph, int km, std::string reads_file, std::string gaf_out, robin_hood::unordered_map<std::string, float>& coverages);
void add_coverages_to_graph(std::string gfa, robin_hood::unordered_map<std::string, float>& coverages);

void compute_exact_CIGARs(std::string gfa_in, std::string gfa_out, int max_overlap, int default_overlap);


namespace std {
    template <>
    class hash<std::pair<int, int>> {
    public:
        size_t operator()(const std::pair<int, bool>& pair) const {
            return std::hash<int>()(pair.first) ^ std::hash<int>()(pair.second);
        }
    };
}


class Segment{
    public:

        std::string name;
        int ID;
        std::vector<std::pair<std::vector<std::pair<int,int>>, std::vector<std::string>>>  links; //first vector is for the links to the left, second vector is for the links to the right. Each link is the ID of the neighbor and its end (0 for left, 1 for right) and the CIGAR

        Segment(){};

        Segment(std::string name, int ID, long int pos_in_file, int length, double coverage){
            this->name = name;
            this->ID = ID;
            this->pos_in_file = pos_in_file;
            this->length = length;
            this->coverage = coverage;
            this->original_coverage = coverage;
            this->haploid = false;
            this->links = std::vector<std::pair<std::vector<std::pair<int,int>>, std::vector<std::string>>>(2);
        }

        Segment(std::string name, int ID, std::vector<std::pair<std::vector<std::pair<int,int>>, std::vector<std::string>>> links, long int pos_in_file, int length, double coverage){
            this->name = name;
            this->ID = ID;
            this->links = links;
            this->pos_in_file = pos_in_file;
            this->length = length;
            this->coverage = coverage;
            this->original_coverage = coverage;
            this->haploid = false;
            this->links = std::vector<std::pair<std::vector<std::pair<int,int>>, std::vector<std::string>>>(2);
        }

        Segment(std::string name, int ID, std::vector<std::pair<std::vector<std::pair<int,int>>, std::vector<std::string>>> links, long int pos_in_file, std::string seq, int length, double coverage){
            this->name = name;
            this->ID = ID;
            this->links = links;
            this->pos_in_file = pos_in_file;
            this->seq = seq;
            this->length = length;
            this->coverage = coverage;
            this->original_coverage = coverage;
            this->haploid = false;
            this->links = std::vector<std::pair<std::vector<std::pair<int,int>>, std::vector<std::string>>>(2);
        }

        bool is_haploid(){return this->haploid;}
        std::vector<std::pair<int,bool>> get_consensus_left(){return this->consensus_left;}
        std::vector<std::pair<int,bool>> get_consensus_right(){return this->consensus_right;}
        long int get_pos_in_file(){return this->pos_in_file;}
        double get_coverage(){return this->coverage;}
        double get_original_coverage(){return this->original_coverage;}
        int get_length(){return this->length;}

        std::string get_seq(std::string& gfa_file){
            
            if (seq != ""){
                // cout << "the seq is already loaded" << endl;
                return seq;
            }

            std::ifstream gfa(gfa_file);
            gfa.seekg(pos_in_file);
            std::string line;
            std::getline(gfa, line);
            std::stringstream ss(line);
            std::string nothing, name, seq;
            ss >> nothing >> name >> seq;
            gfa.close();
            if (seq[0] != 'D'){ //in case of contigs of length 0, we dont want to return "DP:f:54.045"
                return seq;
            }
            else{
                return "";
            }
        }

        void add_neighbor(std::vector<std::pair<int,bool>> new_neighbor, bool left);

        void decrease_coverage(double coverage_out){
            coverage -= coverage_out;
            if (coverage < 1){
                coverage = 1;
            }
        }

        void compute_consensuses();

        std::vector<std::vector<std::pair<int,bool>>> get_strong_neighbors_left(int min_coverage);
        std::vector<std::vector<std::pair<int,bool>>> get_strong_neighbors_right(int min_coverage);

        std::vector<std::vector<std::pair<int,bool>>> get_neighbors_left(){return neighbors_left;}
        std::vector<std::vector<std::pair<int,bool>>> get_neighbors_right(){return neighbors_right;}

        std::vector<std::vector<std::pair<std::pair<int,bool>,std::pair<int,int>>> > get_neighbors_left_with_strengths(){return neighbors_left_with_strengths;}
        std::vector<std::vector<std::pair<std::pair<int,bool>,std::pair<int,int>>> > get_neighbors_right_with_strengths(){return neighbors_right_with_strengths;}

        bool operator!=(const Segment& other) const {
            return !(*this == other);
        }

        bool operator==(const Segment& other) const {
            return name == other.name &&
                   ID == other.ID &&
                   links == other.links &&
                   pos_in_file == other.pos_in_file &&
                   seq == other.seq &&
                   length == other.length &&
                   coverage == other.coverage &&
                   original_coverage == other.original_coverage &&
                   haploid == other.haploid &&
                   consensus_left == other.consensus_left &&
                   consensus_right == other.consensus_right &&
                   neighbors_left_with_strengths == other.neighbors_left_with_strengths &&
                   neighbors_right_with_strengths == other.neighbors_right_with_strengths &&
                   neighbors_left == other.neighbors_left &&
                   neighbors_right == other.neighbors_right;
        }

        //the hash of the segment is the hash of the name
        size_t hash() const{
            return std::hash<std::string>{}(name);
        }

        std::string seq;


    private:
        //consensus sequences of contigs left and right
        std::vector<std::pair<int,bool>> consensus_left;
        std::vector<std::pair<int,bool>> consensus_right;

        std::vector<std::vector<std::pair<std::pair<int,bool>,std::pair<int,int>>>> neighbors_left_with_strengths; //first pair is the path, second pair is the number of reads supporting and disagreeing with the path
        std::vector<std::vector<std::pair<std::pair<int,bool>,std::pair<int,int>>>> neighbors_right_with_strengths; //first pair is the path, second pair is the number of reads supporting and disagreeing with the path

        std::vector<std::vector<std::pair<int,bool>>> neighbors_left; //each path is a std::vector of std::pairs (ID, orientation) of the contigs in the path and a strength (number of reads)
        std::vector<std::vector<std::pair<int,bool>>> neighbors_right; //each path is a std::vector of std::pairs (ID, orientation) of the contigs in the path and a strength (number of reads)

        long int pos_in_file;
        double coverage;
        double original_coverage; //same thing as coverage but cannot be decreased
        int length;
        bool haploid;

};

void load_GFA(std::string gfa_file, std::vector<Segment> &segments, robin_hood::unordered_map<std::string, int> &segment_IDs, bool load_in_RAM);
void merge_adjacent_contigs(std::vector<Segment> &old_segments, std::vector<Segment> &new_segments, std::string original_gfa_file, bool rename, int num_threads);
void output_graph(std::string gfa_output, std::string gfa_input, std::vector<Segment> &segments);



#endif