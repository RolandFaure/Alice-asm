#ifndef REDUCE_AND_EXPAND_H
#define REDUCE_AND_EXPAND_H


#include <string>
#include <unordered_map>
#include <set>
#include "robin_hood.h"

void split_reads_on_non_ACGT_characters(std::string original_read_file, std::string new_read_file);

void reduce(std::string input_file, std::string output_file, int context_length, int compression, int num_threads, bool homopolymer_compression);

void list_kmers_needed_for_expansion(std::string asm_reduced, int km, std::set<std::string> &kmers_needed);
void go_through_the_reads_again_and_index_interesting_kmers(std::string reads_file, std::string assemblyFile, int context_length, int compression, int km, std::set<std::string> &kmers_in_assembly, robin_hood::unordered_map<std::string, std::pair<unsigned long long,unsigned long long>>& kmers, std::string kmer_file, int num_threads, bool homopolymer_compression);
void expand(std::string asm_reduced, std::string output, int km, std::string kmers_file, robin_hood::unordered_map<std::string, std::pair<unsigned long long,unsigned long long>>& kmers);

#endif