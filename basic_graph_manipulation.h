#ifndef BGM_H
#define BGM_H

#include <string>
#include <unordered_map>
#include "robin_hood.h"
#include "rolling_hash.h"

std::string reverse_complement(std::string& seq);
void gfa_to_fasta(std::string gfa, std::string fasta);
void sort_GFA(std::string gfa);

void shave(std::string input_file, std::string output_file, int max_length);
void pop_and_shave_graph(std::string gfa_in, int abundance_min, int min_length, bool contiguity, int k, std::string gfa_out);
void merge_adjacent_contigs_BCALM(std::string gfa_in, std::string gfa_out, int k, std::string path_to_bcalm, std::string path_convertToGFA, std::string path_tmp_folder);

void create_gaf_from_unitig_graph(std::string unitig_graph, int km, std::string reads_file, std::string gaf_out, robin_hood::unordered_map<std::string, float>& coverages);
void add_coverages_to_graph(std::string gfa, robin_hood::unordered_map<std::string, float>& coverages);

void compute_exact_CIGARs(std::string gfa_in, std::string gfa_out, int max_overlap, int default_overlap);

#endif