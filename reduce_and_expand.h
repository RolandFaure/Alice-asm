#ifndef REDUCE_AND_EXPAND_H
#define REDUCE_AND_EXPAND_H


#include <string>
#include <unordered_map>
#include "robin_hood.h"

void reduce(std::string input_file, std::string output_file, int context_length, int compression, int num_threads);
void go_through_the_reads_again(std::string reads_file, std::string assemblyFile, int context_length, int compression, int km, robin_hood::unordered_map<std::string, std::pair<std::string,std::string>>& kmers, int num_threads);
void expand(std::string asm_reduced, std::string output, int km, int length_of_overlaps, robin_hood::unordered_map<std::string, std::pair<std::string,std::string>>& kmers);

#endif