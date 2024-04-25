#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <string>
#include <iostream>

void assembly_hifiasm(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_hifiasm, std::string parameters);
void assembly_bcalm(std::string read_file, int min_abundance, bool contiguity, std::string tmp_folder, int num_threads, std::string final_gfa, std::string path_to_bcalm, std::string path_convertToGFA, std::string path_graphunzip, std::string parameters);
void assembly_spades(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_spades, std::string parameters);
void assembly_minia(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_gatb, std::string path_convertToGFA, std::string parameters);
void assembly_raven(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_raven, std::string parameters);
void assembly_flye(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_flye, std::string parameters);
void assembly_miniasm(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_miniasm, std::string path_to_minimap2, std::string path_to_minipolish, std::string parameters);
void assembly_megahit(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_megahit, std::string path_fastg2gfa, std::string parameters);

#endif