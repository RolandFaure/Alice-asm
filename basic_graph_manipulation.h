#ifndef BGM_H
#define BGM_H

#include <string>
#include <unordered_map>
#include "robin_hood.h"
#include "rolling_hash.h"

std::string reverse_complement(std::string& seq);
void shave(std::string input_file, std::string output_file, int max_length);
void pop_and_shave_homopolymer_errors(std::string gfa_in, std::string gfa_out);
void compute_exact_CIGARs(std::string gfa_in, std::string gfa_out, int max_overlap);
void sort_GFA(std::string gfa);
void create_gaf_from_unitig_graph(std::string unitig_graph, int km, std::string reads_file, std::string gaf_out, robin_hood::unordered_map<std::string, float>& coverages);
void create_gaf_from_unitig_graph2(std::string unitig_graph, int km, std::string reads_file, std::string gaf_out, robin_hood::unordered_map<std::string, float>& coverages);
void merge_adjacent_contigs_BCALM(std::string gfa_in, std::string gfa_out, int k, std::string path_to_bcalm);

void add_coverages_to_graph(std::string gfa, robin_hood::unordered_map<std::string, float>& coverages);

void gfa_to_fasta(std::string gfa, std::string fasta);

#endif