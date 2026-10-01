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
    string command_hifiasm = path_to_hifiasm + " -o " + shell_quote(hifiasm_output + "hifiasm") + " -t " + std::to_string(num_threads) + " " + shell_quote(read_file) + " " + parameters + " > " + shell_quote(tmp_folder + "hifiasm.log") + " 2>&1";
    string untangled_gfa = hifiasm_output + "hifiasm.p_ctg.gfa";

    auto hifiasm_ok = system(command_hifiasm.c_str());
    if (hifiasm_ok != 0){
        cerr << "ERROR: hifiasm failed after running command line\n";
        cerr << command_hifiasm << endl;
        exit(1);
    }

    //move the output to the final file
    std::filesystem::rename(untangled_gfa, final_file);
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
* brief Corrrect reads by building a DBG, trimming and popping bubbles & realigning reads on it
 */
void correct_reads(std::string read_file, int min_abundance, std::string tmp_folder, int num_threads, std::string corrected_reads_file, bool single_genome, std::string path_to_bcalm, std::string path_convertToGFA){
    
    time_t now2 = time(0);
    tm *ltm2 = localtime(&now2);

    cout << " - Correcting reads [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    int kmer_len = 25;

    // launch bcalm        
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << "    - Unitig generation with bcalm [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    int abundancemin = 2;
    if (single_genome){ //if you have a single genome, aggressively delete low coverage kmers
        abundancemin = min_abundance;
    }
    string assembly_file = tmp_folder+"bcalm_correction"+std::to_string(kmer_len);
    string bcalm_command = path_to_bcalm + " -in " + shell_quote(read_file) + " -kmer-size "+std::to_string(kmer_len)+" -abundance-min "+std::to_string(abundancemin)+" -nb-cores "+std::to_string(num_threads)
        + " -out "+shell_quote(assembly_file)+" > "+shell_quote(tmp_folder+"bcalm.log")+" 2>&1";
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
    cout << "    - Converting result to GFA [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    string unitig_file_fa = assembly_file + ".unitigs.fa";
    string unitig_file_gfa = assembly_file + ".unitigs.gfa";
    string convert_command = path_convertToGFA + " " + shell_quote(unitig_file_fa) + " " + shell_quote(unitig_file_gfa) +" "+ std::to_string(kmer_len) + " > " + shell_quote(tmp_folder + "convertToGFA.log") + " 2>&1";
    if (system(convert_command.c_str()) != 0){
        cerr << "ERROR: convertToGFA failed\n";
        cout << convert_command << endl;
        exit(1);
    }
    auto time_convert = std::chrono::high_resolution_clock::now();

    // shave the resulting graph
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << "    - Shaving the graph of small dead ends [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    string shaved_gfa = assembly_file+".unitigs.shaved.gfa";
    pop_and_shave_graph(unitig_file_gfa, min_abundance, 2*kmer_len+10, kmer_len, shaved_gfa, 0, num_threads, single_genome);
    auto time_shave = std::chrono::high_resolution_clock::now();

    //merge the adjacent contigs
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << "    - Merging resulting contigs [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    string merged_gfa = assembly_file+".unitigs.shaved.merged.gfa";
    unordered_map<string, int> segments_IDs2;
    vector<Segment> segments2;
    vector<Segment> merged_segments2;
    load_GFA(shaved_gfa, segments2, segments_IDs2, true); //the last true to load the contigs in memory
    merge_adjacent_contigs(segments2, merged_segments2, shaved_gfa, true, num_threads); //the last bool is to rename the contigs
    output_graph(merged_gfa, shaved_gfa, merged_segments2);
    auto time_merge = std::chrono::high_resolution_clock::now();

    //sort the gfa to have S lines before L lines
    now2 = time(0);
    ltm2 = localtime(&now2);
    sort_GFA(merged_gfa);

    auto time_sort = std::chrono::high_resolution_clock::now();

    //untangle the graph to improve contiguity
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << "    - Aligning the reads to the graph [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    robin_hood::unordered_flat_map<string,float> coverages;
    bool hard_correct = true;  //if hard_correct is true, only keep the reads that can be perfectly corrected (i.e. for which we can find a single path in the graph). If false, keep all reads but only correct the part of the read that can be unambiguously corrected
    create_corrected_reads_from_unitig_graph(merged_gfa, kmer_len, read_file, corrected_reads_file, hard_correct, coverages, num_threads);   
    auto time_gaf = std::chrono::high_resolution_clock::now(); 
    now2 = time(0);
    ltm2 = localtime(&now2);

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
 * @param path_graphunzip Path to the graphunzip executable
 * @return the value of k of the graph that was kept (>= 31)
 */
int assembly_custom(std::string read_file, int min_abundance, std::string tmp_folder, int num_threads, std::string final_gfa, std::vector<int> kmer_sizes_vector, bool single_genome, std::string path_to_bcalm, std::string path_convertToGFA, std::string path_graphunzip, bool hard_contiguity){
    
    string corrected_reads_file = tmp_folder + "corrected_reads.fa";
    correct_reads(read_file, min_abundance, tmp_folder, num_threads, corrected_reads_file, single_genome, path_to_bcalm, path_convertToGFA);
    read_file = corrected_reads_file;

    time_t now2 = time(0);
    tm *ltm2 = localtime(&now2);

    cout << " - Iterative DBG assemby of the compressed reads with increasing k [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;

    vector<int> values_of_k = kmer_sizes_vector; //size of the kmer used to build the graph (min >= km)
    int round = 0; 
    //input of bcalm: the reads, plus the unitigs of the previous k (bcalm takes a comma-separated list of files, no need to concatenate them)
    string bcalm_input = read_file;
    string unitig_file_gfa, unitig_file_fa, merged_gfa;
    for (auto kmer_len: values_of_k){
        // launch bcalm        
        cout << "    - Launching assembly with k=" << kmer_len << endl;
        now2 = time(0);
        ltm2 = localtime(&now2);
        cout << "       - Unitig generation with bcalm [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
        int abundancemin = 1;
        if (single_genome){ //if you have a single genome, aggressively delete low coverage kmers
            abundancemin = min_abundance;
        }
        string bcalm_command = path_to_bcalm + " -in " + shell_quote(bcalm_input) + " -kmer-size "+std::to_string(kmer_len)+" -abundance-min "+std::to_string(abundancemin)+" -nb-cores "+std::to_string(num_threads)
            + " -out "+shell_quote(tmp_folder+"bcalm"+std::to_string(kmer_len))+" > "+shell_quote(tmp_folder+"bcalm"+std::to_string(kmer_len)+".log")+" 2>&1";
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
        cout << "       - Converting result to GFA [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
        unitig_file_fa = tmp_folder+"bcalm"+std::to_string(kmer_len)+".unitigs.fa";
        unitig_file_gfa = tmp_folder+"bcalm"+std::to_string(kmer_len)+".unitigs.gfa";
        string convert_command = path_convertToGFA + " " + shell_quote(unitig_file_fa) + " " + shell_quote(unitig_file_gfa) +" "+ std::to_string(kmer_len) + " > " + shell_quote(tmp_folder + "convertToGFA.log") + " 2>&1";
        if (system(convert_command.c_str()) != 0){
            cerr << "ERROR: convertToGFA failed\n";
            cout << convert_command << endl;
            exit(1);
        }
        auto time_convert = std::chrono::high_resolution_clock::now();

        merged_gfa = unitig_file_gfa;


        //take the unitigs and put them twice in a fasta file, to be used along with the reads to relaunch the assembly with the next k
        if (round < values_of_k.size()-1){
            cout << "       - Concatenating the contigs to the reads to relaunch assembly with higher k [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
        
            string file_with_higher_kmers = read_file + ".higher_k.fa";
            output_unitigs_for_next_k(merged_gfa, file_with_higher_kmers, values_of_k[round+1], 2, num_threads);
            bcalm_input = read_file + "," + file_with_higher_kmers;
        }

        auto time_nextk = std::chrono::high_resolution_clock::now();

        cout << "       - Times: bcalm " << std::chrono::duration_cast<std::chrono::seconds>(time_bcalm - time_start).count() << "s, convert " << std::chrono::duration_cast<std::chrono::seconds>(time_convert - time_bcalm).count() 
            << "s, output for next k: " << std::chrono::duration_cast<std::chrono::seconds>(time_nextk - time_convert).count() << "s"<< endl;

        round++;
    }

    //among all the assemblies with different k >= 31 (needed for the expansion), keep ideally the one with highest k, but
    //not if its graph is much smaller than the largest one (typically if coverage is too low to build a good graph with the higher k)
    string best_gfa = "";
    long int largest_size = 0;
    int best_k = -1;
    for (auto kmer_len: values_of_k){
        if (kmer_len >= 31){
            string gfa_file = tmp_folder+"bcalm"+std::to_string(kmer_len)+".unitigs.gfa";
            largest_size = std::max(largest_size, (long int) std::filesystem::file_size(gfa_file));
        }
    }
    for (auto kmer_len: values_of_k){
        if (kmer_len >= 31){
            string gfa_file = tmp_folder+"bcalm"+std::to_string(kmer_len)+".unitigs.gfa";
            if (std::filesystem::file_size(gfa_file) > 0.9*largest_size){
                best_gfa = gfa_file;
                best_k = kmer_len;
            }
        }
    }
    if (best_k == -1){
        cerr << "ERROR: no k >= 31 in the list of k values\n";
        exit(1);
    }
    merged_gfa = best_gfa;
    cout << " - Best kmer size is " << best_k;
    if (best_k != values_of_k[values_of_k.size()-1]){
        cout << " above this k the assembly completeness decreases";
    }
    cout << endl;
    cout << " =>Done with the iterative assembly, the graph is in " << merged_gfa << "\n" << endl;

    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << " - Untangling the final compressed assembly [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;

    cout << "    - Improving contiguity of assembly by keeping only most covered paths [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    auto time_start = std::chrono::high_resolution_clock::now();
    string shaved_and_popped_gfa = tmp_folder+"bcalm.unitigs.shaved.popped.gfa";
    cut_links_for_contiguity(merged_gfa, shaved_and_popped_gfa, best_k);
    unordered_map<string, int> segments_IDs;
    vector<Segment> segments;
    vector<Segment> merged_segments;
    load_GFA(shaved_and_popped_gfa, segments, segments_IDs, true);  //the last true to load the contigs in memory
    string shaved_and_popped_merged = tmp_folder+"bcalm.unitigs.shaved.popped.merged.gfa";
    merge_adjacent_contigs(segments, merged_segments, shaved_and_popped_gfa, true, num_threads); //the last bool is to rename the contigs
    output_graph(shaved_and_popped_merged, shaved_and_popped_gfa, merged_segments);
    merged_gfa = shaved_and_popped_merged;
    auto time_pop = std::chrono::high_resolution_clock::now();

    //sort the gfa to have S lines before L lines
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << "    - Sorting the GFA [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    sort_GFA(merged_gfa);

    auto time_sort = std::chrono::high_resolution_clock::now();

    //untangle the graph to improve contiguity
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << "    - Aligning the reads to the graph [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    string gaf_file = tmp_folder+"bcalm.unitigs.shaved.merged.unzipped.gaf";
    robin_hood::unordered_flat_map<string,float> coverages;
    create_gaf_from_unitig_graph(merged_gfa, best_k, read_file, gaf_file, coverages, num_threads);   
    auto time_gaf = std::chrono::high_resolution_clock::now(); 
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << "    - Untangling the graph with GraphUnzip [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    
    string unzipped_gfa = tmp_folder+"bcalm.unitigs.shaved.merged.unzipped.gfa";
    string command_unzip = path_graphunzip + " " + shell_quote(merged_gfa) + " " + shell_quote(gaf_file) + " 5 " + 
        std::to_string(num_threads) + " 0 " + shell_quote(unzipped_gfa) + " 1 " + std::to_string(single_genome)
         + " " + std::to_string(best_k) + " " + shell_quote(tmp_folder + "graphunzip.log");
    cout << "    - Command of graphunzip : " << command_unzip << endl;
    auto unzip_ok = system(command_unzip.c_str());
    if (unzip_ok != 0){
        cerr << "ERROR: unzip failed\n";
        exit(1);
    }
    auto time_unzip = std::chrono::high_resolution_clock::now();

    //trim the tips and isolated contigs that result from the unzipping of the graph. Then merge the adjacent contigs
    string tmp_gfa = tmp_folder+"tmp.gfa";
    trim_tips_isolated_contigs_and_bubbles(unzipped_gfa, min_abundance, 2*best_k, tmp_gfa, single_genome, hard_contiguity);
    segments_IDs.clear();
    segments.clear();
    merged_segments.clear();
    load_GFA(tmp_gfa, segments, segments_IDs, true);
    merge_adjacent_contigs(segments, merged_segments, tmp_gfa, true, num_threads); //the last bool is to rename the contigs
    output_graph(final_gfa, tmp_gfa, merged_segments);
    auto time_trim = std::chrono::high_resolution_clock::now();


    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << " => Done untangling the graph, the final compressed graph is in " << final_gfa << " [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]\n" << endl;

    return best_k;
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
    string command_spades = path_to_spades + " -o " + shell_quote(spades_output + "spades") + " --only-assembler -t " + std::to_string(num_threads) + " -s " + shell_quote(read_file) + " " + parameters + " > " + shell_quote(tmp_folder + "spades.log") + " 2>&1";
    string spades_gfa = spades_output + "spades/assembly_graph_with_scaffolds.gfa";

    auto spades_ok = system(command_spades.c_str());
    if (spades_ok != 0){
        cerr << "ERROR: spades failed after running command line\n";
        cerr << command_spades << endl;
        exit(1);
    }

    //copy the output to the final file
    std::filesystem::copy_file(spades_gfa, final_file, std::filesystem::copy_options::overwrite_existing);
}

void assembly_minia(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_gatb, std::string path_convertToGFA, std::string parameters){

    //recover the absolute path to the tmp_folder
    tmp_folder = std::filesystem::absolute(tmp_folder).string();
    read_file = std::filesystem::absolute(read_file).string();
    final_file = std::filesystem::absolute(final_file).string();

    //rm everything starting with minia in the tmp_folder
    for (const auto& entry : std::filesystem::directory_iterator(tmp_folder)){
        if (entry.path().filename().string().rfind("minia", 0) == 0){
            std::filesystem::remove_all(entry.path());
        }
    }

    string minia_output = tmp_folder + "minia";
    string command_minia = path_gatb + " --no-scaffolding --no-error-correction -s " + shell_quote(read_file) + " --nb-cores " + std::to_string(num_threads) 
        + " -o " + shell_quote(minia_output) + " " + parameters + " > " + shell_quote(tmp_folder + "minia.log") + " 2>&1";

    auto minia_ok = system(command_minia.c_str());
    if (minia_ok != 0){
        cerr << "ERROR: minia failed after running command line\n";
        cerr << command_minia << endl;
        exit(1);
    }

    string minia_fasta = minia_output + "_final.contigs.fa";
    //convert the fasta to gfa
    string minia_gfa = tmp_folder + "minia.gfa";
    string command_convert = path_convertToGFA + " " + shell_quote(minia_fasta) + " " + shell_quote(minia_gfa) + " 241 > " + shell_quote(tmp_folder + "convertToGFA.log") + " 2>&1";
    if (system(command_convert.c_str()) != 0){
        cerr << "ERROR: convertToGFA failed after running command line\n";
        cerr << command_convert << endl;
        exit(1);
    }

    //copy the output to the final file
    std::filesystem::copy_file(minia_gfa, final_file, std::filesystem::copy_options::overwrite_existing);
}

void assembly_raven(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_raven, std::string parameters){
    
    string command_raven = path_to_raven + " --graphical-fragment-assembly " + shell_quote(final_file) + " -t " + std::to_string(num_threads) + " " + shell_quote(read_file) + " " + parameters + " > " + shell_quote(tmp_folder + "raven.log") + " 2>&1";

    auto raven_ok = system(command_raven.c_str());
    if (raven_ok != 0){
        cerr << "ERROR: raven failed after running command line\n";
        cerr << command_raven << endl;
        exit(1);
    }
}

void assembly_megahit(std::string read_file, std::string tmp_folder, int num_threads, std::string final_file, std::string path_to_megahit, std::string path_fastg2gfa, std::string parameters){
    
    //remove a potential already existing megahit folder
    std::filesystem::remove_all(tmp_folder + "megahit");


    string command_megahit = path_to_megahit + " -t " + std::to_string(num_threads) + " -o " + shell_quote(tmp_folder + "megahit") + " -r " + shell_quote(read_file) + " " + parameters + " > " + shell_quote(tmp_folder + "megahit.log") + " 2>&1";
    auto megahit_ok = system(command_megahit.c_str());
    cout << "command_megahit: " << command_megahit << "\n";
    if (megahit_ok != 0){
        cerr << "ERROR: megahit failed after running command line\n";
        cerr << command_megahit << endl;
        exit(1);
    }

    //convert the last intermediate assembly (k141) to fastg then to gfa
    string command_to_fastg = "megahit_toolkit contig2fastg 141 " + shell_quote(tmp_folder + "megahit/intermediate_contigs/k141.contigs.fa") + " > " + shell_quote(tmp_folder + "megahit/intermediate_contigs/k141.contigs.fastg");
    cout << "command_to_fastg: " << command_to_fastg << "\n";
    auto to_fastg_ok = system(command_to_fastg.c_str());
    if (to_fastg_ok != 0){
        cerr << "ERROR: megahit_toolkit contig2fastg failed after running command line\n";
        cerr << command_to_fastg << endl;
        exit(1);
    }

    string command_to_gfa = path_fastg2gfa + " " + shell_quote(tmp_folder + "megahit/intermediate_contigs/k141.contigs.fastg") + " > " + shell_quote(final_file);
    auto to_gfa_ok = system(command_to_gfa.c_str());
    cout << "command_to_gfa: " << command_to_gfa << "\n";
    if (to_gfa_ok != 0){
        cerr << "ERROR: fastg2gfa failed after running command line\n";
        cerr << command_to_gfa << endl;
        exit(1);
    }
}

