#include "assembly.h"
#include "basic_graph_manipulation.h"
#include "robin_hood.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

using std::string;
using std::cerr;
using std::endl;
using std::cout;
using std::vector;
using std::ofstream;
using std::ifstream;
using robin_hood::unordered_map;

/**
 * @brief Assemble the read file with hifiasm and output the final assembly in final_file
 * 
 * @param read_file 
 * @param tmp_folder 
 * @param num_threads 
 * @param final_file 
 */
void assembly_hifiasm(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_hifiasm){
    string hifiasm_output = tmp_folder;
    string command_hifiasm = path_to_hifiasm + " -o " + hifiasm_output + "hifiasm -t " + std::to_string(num_threads) + " " + read_file + " > " + tmp_folder + "hifiasm.log 2>&1";
    string untangled_gfa = hifiasm_output + "hifiasm.bp.p_ctg.gfa";

    auto hifiasm_ok = system(command_hifiasm.c_str());
    if (hifiasm_ok != 0){
        cerr << "ERROR: hifiasm failed after running command line\n";
        cerr << command_hifiasm << endl;
        exit(1);
    }

    //move the output to the final file
    string command_move = "mv " + untangled_gfa + " " + final_file;
    system(command_move.c_str());
}

/**
 * @brief Assemble the read file with k-iterative bcalm + graphunzip and output the final assembly in final_file
 * 
 * @param read_file Input read file
 * @param min_abundance Minimum abundance of kmers to consider a kmers as valid
 * @param tmp_folder Folder to store temporary files
 * @param num_threads Number of threads to use
 * @param final_gfa Output final assembly
 * @param path_to_bcalm Path to the bcalm executable
 * @param path_convertToGFA Path to the convertToGFA executable
 * @param path_src Path to the src folder (to get GraphUnzip)
 */
void assembly_bcalm(std::string read_file, int min_abundance, std::string tmp_folder, int num_threads, std::string final_gfa, std::string path_to_bcalm, std::string path_convertToGFA, std::string path_graphunzip){

    cout << " - Iterative DBG assemby of the compressed reads with increasing k\n";

    string merged_gfa = tmp_folder+"bcalm.unitigs.shaved.merged.gfa";
    vector<int> values_of_k = {16,31,41,71}; //size of the kmer used to build the graph (min >= km)
    for (auto kmer_len: values_of_k){
        // launch bcalm        
        cout << "    - Launching assembly with k=" << kmer_len << endl;
        cout << "       - Unitig generation with bcalm" << endl;
        string bcalm_command = path_to_bcalm + " -in " + read_file + " -kmer-size "+std::to_string(kmer_len)+" -abundance-min " 
            + std::to_string(min_abundance) + " -out "+tmp_folder+"bcalm > "+tmp_folder+"bcalm.log 2>&1";
        auto bcalm_ok = system(bcalm_command.c_str());
        if (bcalm_ok != 0){
            cerr << "ERROR: bcalm failed\n";
            cout << bcalm_command << endl;
            exit(1);
        }

        // convert to gfa
        cout << "       - Converting result to GFA" << endl;
        string unitig_file_fa = tmp_folder+"bcalm.unitigs.fa";
        string unitig_file_gfa = tmp_folder+"bcalm.unitigs.gfa";
        string convert_command = path_convertToGFA + " " + unitig_file_fa + " " + unitig_file_gfa +" "+ std::to_string(kmer_len) + " > " + tmp_folder + "convertToGFA.log 2>&1";
        system(convert_command.c_str());

        // shave the resulting graph
        cout << "       - Shaving the graph of small dead ends" << endl;
        string shaved_gfa = tmp_folder+"bcalm.unitigs.shaved.gfa";
        shave(unitig_file_gfa, shaved_gfa, 2*kmer_len-1);

        //merge the adjacent contigs
        cout << "       - Merging resulting contigs" << endl;
        merge_adjacent_contigs_BCALM(shaved_gfa, merged_gfa, kmer_len, path_to_bcalm, path_convertToGFA, tmp_folder);

        //take the contigs of bcalm.unitigs.shaved.merged.unzipped.gfa and put them in a fasta file min_abundance times, and concatenate with compressed_file
        cout << "       - Concatenating the contigs to the reads to relaunch assembly with higher k" << endl;
        
        //open both compressed_file and bcalm.unitigs.shaved.merged.unzipped.gfa
        ofstream input_compressed(read_file, std::ios_base::app);
        ifstream input_graph(merged_gfa);
        string line;
        while (std::getline(input_graph, line))
        {
            if (line[0] == 'S')
            {
                string name;
                string dont_care;
                string sequence;
                std::stringstream ss(line);
                ss >> dont_care >> name >> sequence;
                for (int i = 0 ; i < min_abundance ; i++){
                    input_compressed << ">" << name << "\n";
                    input_compressed << sequence << "\n";
                }
            }
        }
        input_compressed.close();
        input_graph.close();
    }

    cout << " =>Done with the iterative assembly, the graph is in " << merged_gfa << "\n" << endl;

    cout << " - Untangling the final compressed assembly\n";

    //sort the gfa to have S lines before L lines
    cout << "    - Sorting the GFA" << endl;
    sort_GFA(merged_gfa);

    //untangle the graph to improve contiguity
    cout << "    - Aligning the reads to the graph" << endl;
    string gaf_file = tmp_folder+"bcalm.unitigs.shaved.merged.unzipped.gaf";
    unordered_map<string,float> coverages;
    create_gaf_from_unitig_graph(merged_gfa, values_of_k[values_of_k.size()-1], read_file, gaf_file, coverages);
    add_coverages_to_graph(merged_gfa, coverages);
    
    cout << "    - Untangling the graph with GraphUnzip" << endl;
    string command_unzip = path_graphunzip + " unzip -R -l " + gaf_file + " -g " + merged_gfa + " -o " + final_gfa + " > " + tmp_folder + "graphunzip.log 2>&1";
    auto unzip_ok = system(command_unzip.c_str());
    if (unzip_ok != 0){
        cerr << "ERROR: unzip failed\n";
        exit(1);
    }
    cout << " => Done untangling the graph, the final compressed graph is in " << final_gfa << "\n" << endl;
}

/**
 * @brief assemble the read file with spades and output the final assembly in final_file
 * 
 * @param read_file 
 * @param tmp_folder 
 * @param num_threads 
 * @param final_file 
 */
void assembly_spades(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_spades){
    string spades_output = tmp_folder;
    string command_spades = path_to_spades + " -o " + spades_output + "spades --only-assembler -t " + std::to_string(num_threads) + " -s " + read_file;// + " > " + tmp_folder + "spades.log 2>&1";
    string spades_gfa = spades_output + "spades/assembly_graph_with_scaffolds.gfa";

    auto spades_ok = system(command_spades.c_str());
    if (spades_ok != 0){
        cerr << "ERROR: spades failed after running command line\n";
        cerr << command_spades << endl;
        exit(1);
    }

    //move the output to the final file
    string command_move = "cp " + spades_gfa + " " + final_file;
    system(command_move.c_str());
}

void assembly_minia(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_gatb, std::string path_convertToGFA){

    //recover the absolute path to the tmp_folder
    tmp_folder = std::filesystem::absolute(tmp_folder).string();
    read_file = std::filesystem::absolute(read_file).string();
    final_file = std::filesystem::absolute(final_file).string();

    //rm everything starting with minia in the tmp_folder
    string command_rm = "rm -rf " + tmp_folder + "minia*";
    system(command_rm.c_str());

    string minia_output = tmp_folder + "minia";
    string command_minia = path_gatb + " --no-scaffolding --no-error-correction -s " + read_file + " --nb-cores " + std::to_string(num_threads) 
        + " -o " + minia_output + " > " + tmp_folder + "minia.log 2>&1";

    auto minia_ok = system(command_minia.c_str());
    if (minia_ok != 0){
        cerr << "ERROR: minia failed after running command line\n";
        cerr << command_minia << endl;
        exit(1);
    }

    string minia_fasta = minia_output + "_final.contigs.fa";
    //convert the fasta to gfa
    string minia_gfa = tmp_folder + "minia.gfa";
    string command_convert = path_convertToGFA + " " + minia_fasta + " " + minia_gfa + " 241 > " + tmp_folder + "convertToGFA.log 2>&1";

    //move the output to the final file
    string command_move = "cp " + minia_gfa + " " + final_file;
    system(command_move.c_str());
}

void assembly_raven(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_raven){
    
    string command_raven = path_to_raven + " --graphical-fragment-assembly " + final_file + " -t " + std::to_string(num_threads) + " " + read_file + " > " + tmp_folder + "raven.log 2>&1";

    auto raven_ok = system(command_raven.c_str());
    if (raven_ok != 0){
        cerr << "ERROR: raven failed after running command line\n";
        cerr << command_raven << endl;
        exit(1);
    }
}

