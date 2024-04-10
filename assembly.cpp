#include "assembly.h"
#include "basic_graph_manipulation.h"
#include "robin_hood.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

using std::string;
using std::cerr;
using std::endl;
using std::cout;
using std::vector;
using std::ofstream;
using std::ifstream;
using robin_hood::unordered_map;

void assembly_hifiasm(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file){
    string hifiasm_output = tmp_folder;
    string command_hifiasm = "hifiasm -o " + hifiasm_output + "hifiasm -t " + std::to_string(num_threads) + " " + read_file + " > " + tmp_folder + "hifiasm.log 2>&1";
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

void assembly_bcalm(std::string read_file, int min_abundance, std::string tmp_folder, int num_threads, std::string final_gfa, std::string path_to_bcalm, std::string path_convertToGFA, std::string path_src){

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
    string command_unzip = "python " + path_src + "/GraphUnzip/graphunzip.py unzip -R -l " + gaf_file + " -g " + merged_gfa + " -o " + final_gfa + " > " + tmp_folder + "graphunzip.log 2>&1";
    auto unzip_ok = system(command_unzip.c_str());
    if (unzip_ok != 0){
        cerr << "ERROR: unzip failed\n";
        exit(1);
    }
    cout << " => Done untangling the graph, the final compressed graph is in " << final_gfa << "\n" << endl;
}
