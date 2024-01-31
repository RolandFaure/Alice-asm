#ifndef BGM_H
#define BGM_H

#include <string>

std::string reverse_complement(std::string& seq);
void shave(std::string input_file, std::string output_file, int max_length);
void pop_and_shave_homopolymer_errors(std::string gfa_in, std::string gfa_out);
void compute_exact_CIGARs(std::string gfa_in, std::string gfa_out, int max_overlap);
void sort_GFA(std::string gfa);
void create_gaf_from_unitig_graph(std::string unitig_graph, int km, std::string reads_file, std::string gaf_out);

#endif