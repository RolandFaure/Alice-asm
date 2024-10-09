#include "assembly.h"
#include "basic_graph_manipulation.h"
#include "robin_hood.h"
#include "graphunzip.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <chrono>

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
void assembly_hifiasm(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_hifiasm, std::string parameters){
    string hifiasm_output = tmp_folder;
    string command_hifiasm = path_to_hifiasm + " -o " + hifiasm_output + "hifiasm -t " + std::to_string(num_threads) + " " + read_file + " " + parameters + " > " + tmp_folder + "hifiasm.log 2>&1";
    string untangled_gfa = hifiasm_output + "hifiasm.p_ctg.gfa";

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
void assembly_bcalm(std::string read_file, int min_abundance, bool contiguity, int size_longest_read, std::string tmp_folder, int num_threads, std::string final_gfa, std::string path_to_bcalm, std::string path_convertToGFA, std::string path_graphunzip, std::string parameters){
    
    cout << " - Iterative DBG assemby of the compressed reads with increasing k\n";

    string merged_gfa = tmp_folder+"bcalm.unitigs.shaved.merged.gfa";
    vector<int> values_of_k = {17,31}; //size of the kmer used to build the graph (min >= km)
    for (auto kmer_len: values_of_k){
        // launch bcalm        
        cout << "    - Launching assembly with k=" << kmer_len << endl;
        cout << "       - Unitig generation with bcalm" << endl;
        string bcalm_command = path_to_bcalm + " -in " + read_file + " -kmer-size "+std::to_string(kmer_len)+" -abundance-min 2 -nb-cores "+std::to_string(num_threads)
            + " -out "+tmp_folder+"bcalm > "+tmp_folder+"bcalm.log 2>&1";
        auto time_start = std::chrono::high_resolution_clock::now();
        auto bcalm_ok = system(bcalm_command.c_str());
        if (bcalm_ok != 0){
            cerr << "ERROR: bcalm failed\n";
            cout << bcalm_command << endl;
            exit(1);
        }
        auto time_bcalm = std::chrono::high_resolution_clock::now();

        // convert to gfa
        cout << "       - Converting result to GFA" << endl;
        string unitig_file_fa = tmp_folder+"bcalm.unitigs.fa";
        string unitig_file_gfa = tmp_folder+"bcalm.unitigs.gfa";
        string convert_command = path_convertToGFA + " " + unitig_file_fa + " " + unitig_file_gfa +" "+ std::to_string(kmer_len) + " > " + tmp_folder + "convertToGFA.log 2>&1";
        system(convert_command.c_str());
        auto time_convert = std::chrono::high_resolution_clock::now();

        // shave the resulting graph
        cout << "       - Shaving the graph of small dead ends" << endl;
        string shaved_gfa = tmp_folder+"bcalm.unitigs.shaved.gfa";
        pop_and_shave_graph(unitig_file_gfa, -1, 5*kmer_len, contiguity, kmer_len, shaved_gfa);
        auto time_shave = std::chrono::high_resolution_clock::now();

        //merge the adjacent contigs
        cout << "       - Merging resulting contigs" << endl;
        unordered_map<string, int> segments_IDs;
        vector<Segment> segments;
        vector<Segment> merged_segments;
        load_GFA(shaved_gfa, segments, segments_IDs);
        merge_adjacent_contigs(segments, merged_segments, shaved_gfa, true); //the last bool is to rename the contigs
        output_graph(merged_gfa, shaved_gfa, merged_segments);
        auto time_merge = std::chrono::high_resolution_clock::now();

        // merge_adjacent_contigs_BCALM(shaved_gfa, merged_gfa, kmer_len, path_to_bcalm, path_convertToGFA, tmp_folder);

        //take the contigs of bcalm.unitigs.shaved.merged.unzipped.gfa and put them in a fasta file min_abundance times, and concatenate with compressed_file
        cout << "       - Concatenating the contigs to the reads to relaunch assembly with higher k" << endl;
        cout << "       - Times: bcalm " << std::chrono::duration_cast<std::chrono::seconds>(time_bcalm - time_start).count() << "s, convert " << std::chrono::duration_cast<std::chrono::seconds>(time_convert - time_bcalm).count() << "s, shave " << std::chrono::duration_cast<std::chrono::seconds>(time_shave - time_convert).count() << "s, merge " << std::chrono::duration_cast<std::chrono::seconds>(time_merge - time_shave).count() << "s" << endl;
        
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

    auto time_start = std::chrono::high_resolution_clock::now();
    if (contiguity){
        string shaved_and_popped_gfa = tmp_folder+"bcalm.unitigs.shaved.popped.gfa";
        pop_bubbles(merged_gfa, size_longest_read, shaved_and_popped_gfa);
        unordered_map<string, int> segments_IDs;
        vector<Segment> segments;
        vector<Segment> merged_segments;
        load_GFA(shaved_and_popped_gfa, segments, segments_IDs);
        string shaved_and_popped_merged = tmp_folder+"bcalm.unitigs.shaved.popped.merged.gfa";
        merge_adjacent_contigs(segments, merged_segments, shaved_and_popped_gfa, true); //the last bool is to rename the contigs
        output_graph(shaved_and_popped_merged, shaved_and_popped_gfa, merged_segments);
        merged_gfa = shaved_and_popped_merged;
    }
    auto time_pop = std::chrono::high_resolution_clock::now();

    //sort the gfa to have S lines before L lines
    cout << "    - Sorting the GFA" << endl;
    sort_GFA(merged_gfa);

    auto time_sort = std::chrono::high_resolution_clock::now();

    //untangle the graph to improve contiguity
    cout << "    - Aligning the reads to the graph" << endl;
    string gaf_file = tmp_folder+"bcalm.unitigs.shaved.merged.unzipped.gaf";
    unordered_map<string,float> coverages;
    create_gaf_from_unitig_graph(merged_gfa, values_of_k[values_of_k.size()-1], read_file, gaf_file, coverages);   
    auto time_gaf = std::chrono::high_resolution_clock::now(); 
    cout << "    - Untangling the graph with GraphUnzip" << endl;
    
    // string command_unzip = path_graphunzip + " unzip -R -e -l " + gaf_file + " -g " + merged_gfa + " -o " + final_gfa + " -t " + std::to_string(num_threads) + " > " + tmp_folder + "graphunzip.log 2>&1";
    string command_unzip = path_graphunzip + " " + merged_gfa + " " + gaf_file + " 5 1 1 " + final_gfa + " " + std::to_string(contiguity) + " " + tmp_folder + "graphunzip.log";
    cout << "command of graphunzip : " << command_unzip << endl;
    auto unzip_ok = system(command_unzip.c_str());
    if (unzip_ok != 0){
        cerr << "ERROR: unzip failed\n";
        exit(1);
    }
    auto time_unzip = std::chrono::high_resolution_clock::now();

    cout << " => Done untangling the graph, the final compressed graph is in " << final_gfa << "\n" << endl;
    cout << " Times: bcalm " << std::chrono::duration_cast<std::chrono::seconds>(time_pop - time_start).count() << "s, pop " << std::chrono::duration_cast<std::chrono::seconds>(time_sort - time_pop).count() << "s, sort " << std::chrono::duration_cast<std::chrono::seconds>(time_gaf - time_sort).count() << "s, gaf " << std::chrono::duration_cast<std::chrono::seconds>(time_unzip - time_gaf).count() << "s, unzip " << std::chrono::duration_cast<std::chrono::seconds>(time_unzip - time_start).count() << "s" << endl;
}

/**
 * @brief assemble the read file with spades and output the final assembly in final_file
 * 
 * @param read_file 
 * @param tmp_folder 
 * @param num_threads 
 * @param final_file 
 */
void assembly_spades(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_spades, std::string parameters){
    string spades_output = tmp_folder;
    string command_spades = path_to_spades + " -o " + spades_output + "spades --only-assembler -t " + std::to_string(num_threads) + " -s " + read_file + " " + parameters + " > " + tmp_folder + "spades.log 2>&1";
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

void assembly_minia(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_gatb, std::string path_convertToGFA, std::string parameters){

    //recover the absolute path to the tmp_folder
    tmp_folder = std::filesystem::absolute(tmp_folder).string();
    read_file = std::filesystem::absolute(read_file).string();
    final_file = std::filesystem::absolute(final_file).string();

    //rm everything starting with minia in the tmp_folder
    string command_rm = "rm -rf " + tmp_folder + "minia*";
    system(command_rm.c_str());

    string minia_output = tmp_folder + "minia";
    string command_minia = path_gatb + " --no-scaffolding --no-error-correction -s " + read_file + " --nb-cores " + std::to_string(num_threads) 
        + " -o " + minia_output + " " + parameters + " > " + tmp_folder + "minia.log 2>&1";

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

void assembly_raven(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_raven, std::string parameters){
    
    string command_raven = path_to_raven + " --graphical-fragment-assembly " + final_file + " -t " + std::to_string(num_threads) + " " + read_file + " " + parameters + " > " + tmp_folder + "raven.log 2>&1";

    auto raven_ok = system(command_raven.c_str());
    if (raven_ok != 0){
        cerr << "ERROR: raven failed after running command line\n";
        cerr << command_raven << endl;
        exit(1);
    }
}

void assembly_flye(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_flye, std::string parameters){
    
    string command_flye = path_to_flye + " --pacbio-raw " + read_file + " --out-dir " + tmp_folder + "flye --threads " + std::to_string(num_threads) + " " + parameters+ " > " + tmp_folder + "flye.log 2>&1";

    auto flye_ok = system(command_flye.c_str());
    if (flye_ok != 0){
        cerr << "ERROR: flye failed after running command line\n";
        cerr << command_flye << endl;
        exit(1);
    }

    //move the output to the final file
    string command_move = "mv " + tmp_folder + "flye/assembly_graph.gfa " + final_file;
    system(command_move.c_str());
}

void assembly_miniasm(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_miniasm, std::string path_to_minimap2, std::string path_to_minipolish, std::string parameters){

    //all-vs-all read mapping
    string command_minimap = path_to_minimap2 + " -t " + std::to_string(num_threads) + " -x ava-ont " + read_file + " " + read_file + " > " + tmp_folder + "minimap.paf 2> " + tmp_folder + "minimap.log";
    auto minimap_ok = system(command_minimap.c_str());
    if (minimap_ok != 0){
        cerr << "ERROR: minimap failed after running command line\n";
        cerr << command_minimap << endl;
        exit(1);
    }

    //miniasm assembly to get the raw, unpolished assembly
    string command_miniasm = path_to_miniasm + " -f " + read_file + " " + tmp_folder + "minimap.paf > " + tmp_folder + "miniasm.gfa 2> " + tmp_folder + "miniasm.log";
    auto miniasm_ok = system(command_miniasm.c_str());
    if (miniasm_ok != 0){
        cerr << "ERROR: miniasm failed after running command line\n";
        cerr << command_miniasm << endl;
        exit(1);
    }

    //minipolish to polish the assembly
    string command_minipolish = path_to_minipolish + " -t " + std::to_string(num_threads) + " " + read_file + " " + tmp_folder + "miniasm.gfa > " + final_file + " 2> " + tmp_folder + "minipolish.log";
    auto minipolish_ok = system(command_minipolish.c_str());
    if (minipolish_ok != 0){
        cerr << "ERROR: minipolish failed after running command line\n";
        cerr << command_minipolish << endl;
        exit(1);
    }


}

void assembly_megahit(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_megahit, std::string path_fastg2gfa, std::string parameters){
    
    //remove a potential already existing megahit folder
    string command_rm = "rm -rf " + tmp_folder + "megahit > /dev/null 2>&1";
    system(command_rm.c_str());


    string command_megahit = path_to_megahit + " -t " + std::to_string(num_threads) + " -o " + tmp_folder + "megahit -r " + read_file + " " + parameters + " > " + tmp_folder + "megahit.log 2>&1";
    auto megahit_ok = system(command_megahit.c_str());
    cout << "command_megahit: " << command_megahit << "\n";
    if (megahit_ok != 0){
        cerr << "ERROR: megahit failed after running command line\n";
        cerr << command_megahit << endl;
        exit(1);
    }

    //convert the last intermediate assembly (k141) to fastg then to gfa
    string command_to_fastg = "megahit_toolkit contig2fastg 141 " + tmp_folder + "megahit/intermediate_contigs/k141.contigs.fa > " + tmp_folder + "megahit/intermediate_contigs/k141.contigs.fastg";
    cout << "command_to_fastg: " << command_to_fastg << "\n";
    auto to_fastg_ok = system(command_to_fastg.c_str());
    if (to_fastg_ok != 0){
        cerr << "ERROR: megahit_toolkit contig2fastg failed after running command line\n";
        cerr << command_to_fastg << endl;
        exit(1);
    }

    string command_to_gfa = path_fastg2gfa + " " + tmp_folder + "megahit/intermediate_contigs/k141.contigs.fastg > " + final_file;
    auto to_gfa_ok = system(command_to_gfa.c_str());
    cout << "command_to_gfa: " << command_to_gfa << "\n";
    if (to_gfa_ok != 0){
        cerr << "ERROR: fastg2gfa failed after running command line\n";
        cerr << command_to_gfa << endl;
        exit(1);
    }
}

