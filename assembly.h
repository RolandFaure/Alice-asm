#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <string>
#include <iostream>

void assembly_hifiasm(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file);
void assembly_bcalm(std::string read_file, int min_abundance, std::string tmp_folder, int num_threads, std::string final_gfa, std::string path_to_bcalm, std::string path_convertToGFA, std::string path_src);

#endif