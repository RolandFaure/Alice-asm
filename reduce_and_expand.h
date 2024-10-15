#ifndef REDUCE_AND_EXPAND_H
#define REDUCE_AND_EXPAND_H


#include <string>
#include <unordered_map>
#include <set>
#include "robin_hood.h"

void reduce(std::string input_file, std::string output_file, int context_length, int compression, int num_threads, bool homopolymer_compression);

void list_kmers_needed_for_expansion(std::string asm_reduced, int km, std::set<uint64_t> &kmers_needed);
void go_through_the_reads_again_and_index_interesting_kmers(std::string reads_file, std::string assemblyFile, int context_length, int compression, int km, std::set<uint64_t> &kmers_in_assembly, robin_hood::unordered_map<uint64_t, std::pair<unsigned long long,unsigned long long>>& kmers, std::string kmer_file, int num_threads, bool homopolymer_compression);
void expand(std::string asm_reduced, std::string output, int km, std::string kmers_file, robin_hood::unordered_map<uint64_t, std::pair<unsigned long long,unsigned long long>>& kmers);

#endif