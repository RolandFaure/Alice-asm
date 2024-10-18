#ifndef REDUCE_AND_EXPAND_H
#define REDUCE_AND_EXPAND_H


#include <string>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include "robin_hood.h"

void reduce(std::string input_file, std::string output_file, int context_length, int compression, int num_threads, bool homopolymer_compression);

void expand_or_list_kmers_needed_for_expansion(std::string mode, std::string asm_reduced, int km, std::unordered_set<uint64_t> &central_kmers_needed, std::unordered_set<uint64_t> &full_kmers_needed, std::string central_kmers_file, std::string full_kmers_file, robin_hood::unordered_map<uint64_t, std::pair<unsigned long long,unsigned long long>>& kmers, std::string output);
void go_through_the_reads_again_and_index_interesting_kmers(std::string reads_file, std::string assemblyFile, int context_length, int compression, int km, std::unordered_set<uint64_t> &central_kmers_in_assembly, std::unordered_set<uint64_t> &full_kmers_in_assembly, robin_hood::unordered_map<uint64_t, std::pair<unsigned long long,unsigned long long>>& kmers, std::string central_kmer_file, std::string full_kmer_file, int num_threads, bool homopolymer_compression);

#endif