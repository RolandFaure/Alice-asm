#ifndef BGM_H
#define BGM_H

#include <string>
#include <array>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include "robin_hood.h"
#include "rolling_hash.h"

std::string reverse_complement(std::string& seq);

//quote a file path for a shell command line (paths may contain spaces or shell characters)
inline std::string shell_quote(const std::string& path){
    std::string quoted = "'";
    for (char c : path){
        if (c == '\''){
            quoted += "'\\''";
        }
        else{
            quoted += c;
        }
    }
    return quoted + "'";
}
void gfa_to_fasta(std::string gfa, std::string fasta);
void sort_GFA(std::string gfa);

void pop_and_shave_graph(std::string gfa_in, int abundance_min, int min_length, int k, std::string gfa_out, int extra_coverage, int num_threads, bool single_genome);
void cut_links_for_contiguity(std::string gfa_in, std::string gfa_out, int k);
void trim_tips_isolated_contigs_and_bubbles(std::string gfa_in, int min_coverage, int min_length, std::string gfa_out, bool single_genome, bool hard_contiguity);

void create_corrected_reads_from_unitig_graph(std::string unitig_graph, int km, std::string reads_file, std::string output_file, bool hard_correct, robin_hood::unordered_flat_map<std::string, float>& coverages, int num_threads);
void create_gaf_from_unitig_graph(std::string unitig_graph, int km, std::string reads_file, std::string output_file, robin_hood::unordered_flat_map<std::string, float>& coverages, int num_threads);
void add_coverages_to_graph(std::string gfa, robin_hood::unordered_map<std::string, float>& coverages);

void compute_exact_CIGARs(std::string gfa_in, std::string gfa_out, int max_overlap, int default_overlap, double bases_per_compressed_base, int num_threads);
long long total_sequence_length(std::string gfa);
int overlap_length_of_CIGAR(const std::string& cigar);


//hash of a pair of ints (e.g. (ID, end) of a segment)
struct PairHash {
    size_t operator()(const std::pair<int, int>& p) const {
        return std::hash<long long>()(((long long) p.first << 32) ^ (unsigned int) p.second);
    }
};


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
            this->links = std::vector<std::pair<std::vector<std::pair<int,int>>, std::vector<std::string>>>(2);
        }

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

        void decrease_coverage(double coverage_out){
            coverage -= coverage_out;
            if (coverage < 1){
                coverage = 1;
            }
        }

        std::string seq;


    private:
        long int pos_in_file;
        double coverage;
        double original_coverage; //same thing as coverage but cannot be decreased
        int length;

};

void load_GFA(std::string gfa_file, std::vector<Segment> &segments, robin_hood::unordered_map<std::string, int> &segment_IDs, bool load_in_RAM);
void merge_adjacent_contigs(std::vector<Segment> &old_segments, std::vector<Segment> &new_segments, std::string original_gfa_file, bool rename, int num_threads);
void output_graph(std::string gfa_output, std::string gfa_input, std::vector<Segment> &segments);





#endif