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
#include <omp.h>

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
 * @brief Output exactly num_copies times all k-mers of the unitigs in the reads_fa file
 * 
 * @param unitig_gfa graph of the unitigs (built with a smaller k)
 * @param reads_fa 
 * @param k 
 * @param num_copies
 * @param bcalm path to the bcalm executable
 * @param num_threads
 */
void output_unitigs_for_next_k(std::string unitig_gfa, std::string file_with_higher_kmers, int k, int num_copies, int num_threads){

    //load all the links of the unitigs
    ifstream gfa(unitig_gfa);
    string line;
    ofstream out(file_with_higher_kmers);
    while (getline(gfa, line)){
        string nothing;
        if (line[0] == 'S'){
            string name;
            string seq;
            std::stringstream ss(line);
            ss >> nothing >> name >> seq;
            if (seq.length() >= k) {
                for (int i = 0; i < num_copies; ++i) {
                    out << ">" << name << "_" << i << "\n" << seq << "\n";
                }
            }
        }
    }
    gfa.close();
    out.close();
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
void assembly_custom(std::string read_file, int min_abundance, bool contiguity, int size_longest_read, std::string tmp_folder, int num_threads, std::string final_gfa, std::vector<int> kmer_sizes_vector, bool single_genome, std::string path_to_bcalm, std::string path_convertToGFA, std::string path_graphunzip){
    
    time_t now2 = time(0);
    tm *ltm2 = localtime(&now2);

    cout << " - Iterative DBG assemby of the compressed reads with increasing k [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;

    vector<int> values_of_k = kmer_sizes_vector; //size of the kmer used to build the graph (min >= km)
    int round = 0; 
    string file_with_unitigs_from_past_k_and_reads = read_file+".with_unitigs_from_previous_k.fa";
    //copy the reads to the file with unitigs
    string command_copy = "cp " + read_file + " " + file_with_unitigs_from_past_k_and_reads;
    system(command_copy.c_str());
    for (auto kmer_len: values_of_k){
        // launch bcalm        
        cout << "    - Launching assembly with k=" << kmer_len << endl;
        now2 = time(0);
        ltm2 = localtime(&now2);
        cout << "       - Unitig generation with bcalm [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
        int abundancemin = 2;
        if (single_genome){ //if you have a single genome, aggressively delete low coverage kmers
            abundancemin = min_abundance;
        }
        string bcalm_command = path_to_bcalm + " -in " + file_with_unitigs_from_past_k_and_reads + " -kmer-size "+std::to_string(kmer_len)+" -abundance-min "+std::to_string(abundancemin)+" -nb-cores "+std::to_string(num_threads)
            + " -out "+tmp_folder+"bcalm"+std::to_string(kmer_len)+" > "+tmp_folder+"bcalm.log 2>&1";
        auto time_start = std::chrono::high_resolution_clock::now();
        auto bcalm_ok = system(bcalm_command.c_str());
        if (bcalm_ok != 0){
            cerr << "ERROR: bcalm failed\n";
            cout << bcalm_command << endl;
            exit(1);
        }
        auto time_bcalm = std::chrono::high_resolution_clock::now();

        // convert to gfa
        now2 = time(0);
        ltm2 = localtime(&now2);
        cout << "       - Converting result to GFA [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
        string unitig_file_fa = tmp_folder+"bcalm"+std::to_string(kmer_len)+".unitigs.fa";
        string unitig_file_gfa = tmp_folder+"bcalm"+std::to_string(kmer_len)+".unitigs.gfa";
        string convert_command = path_convertToGFA + " " + unitig_file_fa + " " + unitig_file_gfa +" "+ std::to_string(kmer_len) + " > " + tmp_folder + "convertToGFA.log 2>&1";
        auto res = system(convert_command.c_str());
        auto time_convert = std::chrono::high_resolution_clock::now();

        // shave the resulting graph and //-min_abundance on all the abundances for each round (to remove the contigs that were added at the end of assembly for higher k)
        now2 = time(0);
        ltm2 = localtime(&now2);
        cout << "       - Shaving the graph of small dead ends [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
        string shaved_gfa = tmp_folder+"bcalm"+std::to_string(kmer_len)+".unitigs.shaved.gfa";
        pop_and_shave_graph(unitig_file_gfa, min_abundance, 2*kmer_len+1, contiguity, kmer_len, shaved_gfa, std::min(1,round)*2, num_threads, single_genome); //std::min(1,round)*2 because we want to remove the contigs that were added at the end of the previous assembly in two copies
        auto time_shave = std::chrono::high_resolution_clock::now();

        //merge the adjacent contigs
        now2 = time(0);
        ltm2 = localtime(&now2);
        cout << "       - Merging resulting contigs [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
        string merged_gfa = tmp_folder+"bcalm"+std::to_string(kmer_len)+".unitigs.shaved.merged.gfa";
        if (true || round == values_of_k.size()-1){ //we are in the last round, do a proper merge that keeps the coverages //BUT gfatools does not seem to work well so do this every time...
            unordered_map<string, int> segments_IDs;
            vector<Segment> segments;
            vector<Segment> merged_segments;
            load_GFA(shaved_gfa, segments, segments_IDs, true); //the last true to load the contigs in memory
            merge_adjacent_contigs(segments, merged_segments, shaved_gfa, true, num_threads); //the last bool is to rename the contigs
            output_graph(merged_gfa, shaved_gfa, merged_segments);
        }
        // else{ //merge using gfatools asm -u: much faster, though it does not keep the coverages
        //     string merged_gfa_tmp = tmp_folder+"bcalm"+std::to_string(kmer_len)+".unitigs.shaved.merged.gfa";
        //     string merge_command = "gfatools asm -u "+shaved_gfa+" > "+merged_gfa + " 2> "+tmp_folder+"gfatools.log";
        //     system(merge_command.c_str());
        // }
        auto time_merge = std::chrono::high_resolution_clock::now();

        //take the contigs of bcalm.unitigs.shaved.merged.unzipped.gfa and put them in a fasta file min_abundance times, and concatenate with compressed_file
        if (round < values_of_k.size()-1){
            cout << "       - Concatenating the contigs to the reads to relaunch assembly with higher k" << endl;
        
            string file_with_higher_kmers = read_file + ".higher_k.fa";
            output_unitigs_for_next_k(merged_gfa, file_with_higher_kmers, values_of_k[round+1], 2, num_threads);
            //concatenate the originial reads with the file_with_higher_kmers to relaunch the assembly
            string command_concatenate = "cat " + read_file + " " + file_with_higher_kmers + " > " + file_with_unitigs_from_past_k_and_reads;
            system(command_concatenate.c_str());
        }

        auto time_nextk = std::chrono::high_resolution_clock::now();

        cout << "       - Times: bcalm " << std::chrono::duration_cast<std::chrono::seconds>(time_bcalm - time_start).count() << "s, convert " << std::chrono::duration_cast<std::chrono::seconds>(time_convert - time_bcalm).count() 
            << "s, shave " << std::chrono::duration_cast<std::chrono::seconds>(time_shave - time_convert).count() << "s, merge " << std::chrono::duration_cast<std::chrono::seconds>(time_merge - time_shave).count()  
            << "s, output for next k: " << std::chrono::duration_cast<std::chrono::seconds>(time_nextk - time_merge).count() << "s"<< endl;

        round++;
    }

    string merged_gfa = tmp_folder+"bcalm"+std::to_string(values_of_k[values_of_k.size()-1])+".unitigs.shaved.merged.gfa";
    cout << " =>Done with the iterative assembly, the graph is in " << merged_gfa << "\n" << endl;

    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << " - Untangling the final compressed assembly [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;

    auto time_start = std::chrono::high_resolution_clock::now();
    if (contiguity){
        string shaved_and_popped_gfa = tmp_folder+"bcalm.unitigs.shaved.popped.gfa";
        pop_bubbles(merged_gfa, size_longest_read, shaved_and_popped_gfa);
        unordered_map<string, int> segments_IDs;
        vector<Segment> segments;
        vector<Segment> merged_segments;
        load_GFA(shaved_and_popped_gfa, segments, segments_IDs, true);  //the last true to load the contigs in memory
        string shaved_and_popped_merged = tmp_folder+"bcalm.unitigs.shaved.popped.merged.gfa";
        merge_adjacent_contigs(segments, merged_segments, shaved_and_popped_gfa, true, num_threads); //the last bool is to rename the contigs
        output_graph(shaved_and_popped_merged, shaved_and_popped_gfa, merged_segments);
        merged_gfa = shaved_and_popped_merged;
    }
    auto time_pop = std::chrono::high_resolution_clock::now();

    //sort the gfa to have S lines before L lines
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << "    - Sorting the GFA [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    sort_GFA(merged_gfa);

    auto time_sort = std::chrono::high_resolution_clock::now();

    //untangle the graph to improve contiguity
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << "    - Aligning the reads to the graph [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    string gaf_file = tmp_folder+"bcalm.unitigs.shaved.merged.unzipped.gaf";
    unordered_map<string,float> coverages;
    create_gaf_from_unitig_graph(merged_gfa, values_of_k[values_of_k.size()-1], read_file, gaf_file, coverages);   
    auto time_gaf = std::chrono::high_resolution_clock::now(); 
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << "    - Untangling the graph with GraphUnzip [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    
    // string command_unzip = path_graphunzip + " unzip -R -e -l " + gaf_file + " -g " + merged_gfa + " -o " + final_gfa + " -t " + std::to_string(num_threads) + " > " + tmp_folder + "graphunzip.log 2>&1";
    string unzipped_gfa = tmp_folder+"bcalm.unitigs.shaved.merged.unzipped.gfa";
    string command_unzip = path_graphunzip + " " + merged_gfa + " " + gaf_file + " 5 " + std::to_string(num_threads) + " 0 " + unzipped_gfa + " " + std::to_string(contiguity) + + " " + std::to_string(single_genome) + " " + tmp_folder + "graphunzip.log";
    cout << "    - Command of graphunzip : " << command_unzip << endl;
    auto unzip_ok = system(command_unzip.c_str());
    if (unzip_ok != 0){
        cerr << "ERROR: unzip failed\n";
        exit(1);
    }
    auto time_unzip = std::chrono::high_resolution_clock::now();

    //trim the tips and isolated contigs that result from the unzipping of the graph. Then merge the adjacent contigs
    string tmp_gfa = tmp_folder+"tmp.gfa";
    trim_tips_isolated_contigs_and_bubbles(unzipped_gfa, min_abundance, 2*values_of_k[values_of_k.size()-1], tmp_gfa);
    unordered_map<string, int> segments_IDs;
    vector<Segment> segments;
    vector<Segment> merged_segments;
    load_GFA(tmp_gfa, segments, segments_IDs, true);
    merge_adjacent_contigs(segments, merged_segments, tmp_gfa, true, num_threads); //the last bool is to rename the contigs
    output_graph(final_gfa, tmp_gfa, merged_segments);
    auto time_trim = std::chrono::high_resolution_clock::now();


    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << " => Done untangling the graph, the final compressed graph is in " << final_gfa << " [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]\n" << endl;
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

