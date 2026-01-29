// author : Roland Faure

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <unordered_set>
#include <chrono>
#include <thread>
#include <omp.h> //for efficient parallelization
#include <set>

#include "reduce_and_expand.h"
#include "basic_graph_manipulation.h"
#include "robin_hood.h"
#include "assembly.h"
#include "clipp.h"
#include "test.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;
using robin_hood::unordered_map;
using std::pair;
using std::cerr;
using std::ifstream;
using std::ofstream;
using std::unordered_set;
using std::set;

// ANSI escape codes for text color
#define RED_TEXT "\033[1;31m"
#define GREEN_TEXT "\033[1;32m"
#define RESET_TEXT "\033[0m"

string version = "0.7.3";
string date = "2026-01-28";
string author = "Roland Faure";

//small function to exaceute a shell command and catch the result
std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    // Remove trailing newline, if present
    if (!result.empty() && result.back() == '\n') {
        result.pop_back();
    }
    return result;
}

void check_dependencies(string assembler, string path_bcalm, string path_hifiasm, string path_spades, string path_minia, string path_raven, string path_to_flye, string path_minimap, string path_miniasm, string path_minipolish, string path_megahit, string path_fastg2gfa,
    string &path_convertToGFA, string &path_graphunzip, string path_src){

        
    string command_bcalm = path_bcalm + " --help 2> trash.log > trash.log";
    auto bcalm_ok = system(command_bcalm.c_str());

    string command_hifiasm = path_hifiasm + " -h 2> trash.log > trash.log";
    auto hifiasm_ok = system(command_hifiasm.c_str());

    string command_spades = path_spades + " --help 2> trash.log > trash.log";
    auto spades_ok = system(command_spades.c_str());

    string command_raven = path_raven + " --help 2> trash.log > trash.log";
    auto raven_ok = system(command_raven.c_str());

    string command_flye = path_to_flye + " --help 2> trash.log > trash.log";
    auto flye_ok = system(command_flye.c_str());

    string command_minimap = path_minimap + " --version 2> trash.log > trash.log";
    auto minimap_ok = system(command_minimap.c_str());

    string command_miniasm = path_miniasm + " -V 2> trash.log > trash.log";
    auto miniasm_ok = system(command_miniasm.c_str());

    string command_minipolish = path_minipolish + " --version 2> trash.log > trash.log";
    auto minipolish_ok = system(command_minipolish.c_str());

    string command_megahit = path_megahit + " --version 2> trash.log > trash.log";
    auto megahit_ok = system(command_megahit.c_str());

    string command_fastg2gfa = path_fastg2gfa + " 2> trash.log > trash.log";
    auto fastg2gfa_ok = system(command_fastg2gfa.c_str());

    string command_python3 = "python3 --version 2> trash.log > trash.log";
    auto python3_ok = system(command_python3.c_str());

    string command_convertToGFA = path_convertToGFA + " -h 2> trash.log > trash.log";
    auto convertToGFA_ok = system(command_convertToGFA.c_str());
    if (convertToGFA_ok != 0) {
        string bad_path = path_convertToGFA;
        path_convertToGFA = "python3 " + exec("which convertToGFA.py");
        convertToGFA_ok = system((path_convertToGFA + " -h 2> trash.log > trash.log").c_str());
        if (convertToGFA_ok != 0 || path_convertToGFA == "python3 ") {
            cerr << "ERROR: convertToGFA.py not found, problem in the installation, error code 321.\n";
            cout << "tried " << endl << bad_path << endl << path_convertToGFA << endl;
            exit(1);
        }
    }


    // string command_minia = path_minia + " --help 2> trash.log > trash.log";
    // auto minia_ok = system(command_minia.c_str());
    // cout << "trying command line " << command_minia << " " << minia_ok << " " << minia_ok2 << endl;

    int graphunzip_ok = 0;
    if (assembler == "custom"){
        auto graphunzip_ok = system((path_graphunzip + " --help >trash.log 2>trash.log ").c_str());
        if (graphunzip_ok != 0){
            path_graphunzip = "graphunzip";
            graphunzip_ok = system((path_graphunzip + " --help >trash.log 2>trash.log ").c_str());
            if (graphunzip_ok != 0){
                cerr << "ERROR: graphunzip not found, problem in the installation, error code 322.\n";
                exit(1);
            }
        }
    }

    // int gfatools_ok = system("gfatools version 2> trash.log > trash.log");

    std::cout << "_______________________________" << std::endl;
    std::cout << "|    Dependency     |  Found  |" << std::endl;
    std::cout << "|-------------------|---------|" << std::endl;
    std::cout << "|    python3        |   " << (python3_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    if (assembler == "custom")
        std::cout << "|    bcalm          |   " << (bcalm_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    else if (assembler == "hifiasm")
        std::cout << "|    hifiasm        |   " << (hifiasm_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    else if (assembler == "spades")
        std::cout << "|    spades         |   " << (spades_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    // else if (assembler == "gatb-minia")
    //     std::cout << "|    gatb-minia     |   " << (minia_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    else if (assembler == "raven")
        std::cout << "|    raven          |   " << (raven_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    else if (assembler == "flye")
        std::cout << "|    flye           |   " << (flye_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    else if (assembler == "miniasm"){
        std::cout << "|    minimap2       |   " << (minimap_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
        std::cout << "|    miniasm        |   " << (miniasm_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
        std::cout << "|    minipolish     |   " << (minipolish_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    }
    else if (assembler == "megahit"){
        std::cout << "|    megahit        |   " << (megahit_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
        std::cout << "|    fastg2gfa      |   " << (fastg2gfa_ok == 256 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    }
    // std::cout << "|    gfatools       |   " << (gfatools_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    std::cout << "-------------------------------" << std::endl;


    if ((bcalm_ok != 0 && assembler == "custom") || 
        (hifiasm_ok != 0 && assembler == "hifiasm") || 
        (spades_ok != 0 && assembler == "spades") || 
        // (minia_ok != 0 && assembler == "gatb-minia") ||
        (raven_ok != 0 && assembler == "raven") ||
        (flye_ok != 0 && assembler == "flye") ||
        ((minimap_ok != 0 || miniasm_ok != 0 || minipolish_ok != 0) && assembler == "miniasm") ||
        (megahit_ok != 0 && assembler == "megahit" && fastg2gfa_ok != 0) ||
        (convertToGFA_ok != 0 || graphunzip_ok != 0) ||
        (python3_ok != 0)){
        std::cout << "Error: some dependencies are missing." << std::endl;
        exit(1);
    }

}

int main(int argc, char** argv)
{
    // unordered_map<string,float> coverages;
    // string a = "out_alice/tmp/bcalm.unitigs.shaved.popped.merged.gfa";
    // string b = "out_alice/tmp/compressed.fa";
    // create_gaf_from_unitig_graph(a, 31, b, "trash.gaf", coverages);
    // exit(0);

    //use clipp to parse the command line
    bool help = false;
    string input_file, output_folder;
    string assembler = "custom";
    string path_to_bcalm = "bcalm";
    string path_to_hifiasm = "hifiasm_meta";
    string path_to_spades = "spades.py";
    string path_to_minia = "gatb";
    string path_to_raven = "raven";
    string path_to_flye = "flye";
    string path_to_miniasm = "miniasm";
    string path_to_minimap2 = "minimap2";
    string path_to_minipolish = "racon";
    string path_to_megahit = "megahit";
    string assembler_parameters = "";
    string test_ref_gfa = "";
    bool rescue = false;
    bool contiguity = true;
    bool single_genome= false;
    int min_abundance = 5;
    string kmer_sizes = "17,21,31,61,101,191";
    int order = 101;
    int compression = 20;
    int num_threads = 1;
    bool no_hpc = false;
    bool clean = false;
    auto cli = (
        //input/output option
        clipp::required("-r", "--reads").doc("input file (fasta/q)") & clipp::opt_value("r",input_file),
        clipp::required("-o", "--output").doc("output folder") & clipp::opt_value("o",output_folder),

        //Performance options
        clipp::option("-t", "--threads").doc("number of threads [1]") & clipp::opt_value("t", num_threads),

        //Compression options
        clipp::option("-l", "--order").doc("order of MSR compression (odd) [101]") & clipp::opt_value("o", order),
        clipp::option("-c", "--compression").doc("compression factor [20]") & clipp::opt_value("c", compression),
        clipp::option("-H", "--no-hpc").set(no_hpc).doc("turn off homopolymer and homodimer compression"),

        //Assembly options for the custom assembler
        clipp::option("-m", "--min-abundance").doc("minimum abundance of kmer to consider solid - RECOMMENDED to set to coverage/2 if single-genome [5]") & clipp::opt_value("m", min_abundance),
        clipp::option("-k", "--kmer-sizes").doc("comma-separated increasing sizes of k for assembly, must go at least to 31 [17,21,31]") & clipp::opt_value("k", kmer_sizes),
        // clipp::option("--contiguity").set(contiguity).doc("Favor contiguity over recovery of rare strains [off]"), //in our tests, contiguity is better turned on
        clipp::option("--single-genome").set(single_genome).doc("Switch on if assembling a single genome"),

        //Other assemblers options
        clipp::option("-a", "--assembler").doc("assembler to use {custom, hifiasm, spades, raven, gatb-minia, megahit} [custom]") & clipp::opt_value("a", assembler),
        
        // clipp::option("--parameters").doc("extra parameters to pass to the assembler (between quotation marks) [\"\"]") & clipp::opt_value("p", assembler_parameters),
        clipp::option("--bcalm").doc("path to bcalm [bcalm]") & clipp::opt_value("b", path_to_bcalm),
        // clipp::option("--hifiasm_meta").doc("path to hifiasm_meta [hifiasm_meta]") & clipp::opt_value("h", path_to_hifiasm),
        clipp::option("--spades").doc("path to spades [spades.py]") & clipp::opt_value("s", path_to_spades),
        // clipp::option("--raven").doc("path to raven [raven]") & clipp::opt_value("r", path_to_bcalm),
        // // clipp::option("--flye").doc("path to flye [flye]") & clipp::opt_value("f", path_to_flye), //flye does not work well with compressed reads
        // clipp::option("--gatb-minia").doc("path to gatb-minia [gatb]") & clipp::opt_value("g", path_to_minia),
        // clipp::option("--megahit").doc("path to megahit [megahit]") & clipp::opt_value("m", path_to_megahit),
        // clipp::option("--miniasm").doc("path to miniasm [miniasm]") & clipp::opt_value("m", path_to_miniasm), //we did not manage to make miniasm work
        // clipp::option("--minimap2").doc("path to minimap2 [minimap2]") & clipp::opt_value("m", path_to_minimap2),
        // clipp::option("--minipolish").doc("path to minipolish [minipolish]") & clipp::opt_value("r", path_to_minipolish),
        
        //Other options
        clipp::option("--clean").set(clean).doc("remove the tmp folder at the end [off]"),
        clipp::option("--test").doc("(developers only) to compare the result against this reference") & clipp::opt_value("t", test_ref_gfa),
        clipp::option("-v", "--version").call([]{ std::cout << "version " << version << "\nLast update: " << date << "\nAuthor: " << author << std::endl; exit(0); }).doc("print version and exit"),
        clipp::option("-h", "--help").set(help).doc("print this help message and exit")
    );

    bool homopolymer_compression = !no_hpc;


    //ascii art of a cake:
    //     _.mnm._    
    //    ( _____ ) 
    //     |     |
    //      `___/    

    //ascii art of a the Alice Assembler:
    //   _______ _                     _ _                                            _     _           
    //  |__   __| |              /\   | (_)              /\                          | |   | |          
    //     | |  | |__   ___     /  \  | |_  ___ ___     /  \   ___ ___  ___ _ __ ___ | |__ | | ___ _ __ 
    //     | |  | '_ \ / _ \   / /\ \ | | |/ __/ _ \   / /\ \ / __/ __|/ _ \ '_ ` _ \| '_ \| |/ _ \ '__|
    //     | |  | | | |  __/  / ____ \| | | (_|  __/  / ____ \\__ \__ \  __/ | | | | | |_) | |  __/ |   
    //     |_|  |_| |_|\___| /_/    \_\_|_|\___\___| /_/    \_\___/___/\___|_| |_| |_|_.__/|_|\___|_|   

    //ascii art of a bottle:
//     ::
//    :  :
//    :  :
//    :__:  


    cout << "   _______ _                     _ _                                            _     _              "           << "      " << "           " << endl;
    cout << "  |__   __| |              /\\   | (_)              /\\                          | |   | |             "         << "      " << "           " << endl;
    cout << "     | |  | |__   ___     /  \\  | |_  ___ ___     /  \\   ___ ___  ___ _ __ ___ | |__ | | ___ _ __    "         << "  ::  " << "   _.mnm._ " << endl;
    cout << "     | |  | '_ \\ / _ \\   / /\\ \\ | | |/ __/ _ \\   / /\\ \\ / __/ __|/ _ \\ '_ ` _ \\| '_ \\| |/ _ \\ '__|   "<< " :  : " << "  ( _____ )" << endl;
    cout << "     | |  | | | |  __/  / ____ \\| | | (_|  __/  / ____ \\\\__ \\__ \\  __/ | | | | | |_) | |  __/ |      "      << " :  : " << "   |     | " << endl;
    cout << "     |_|  |_| |_|\\___| /_/    \\_\\_|_|\\___\\___| /_/    \\_\\___/___/\\___|_| |_| |_|_.__/|_|\\___|_|      "  << " :__: " << "    `___/  " << endl;
    cout << endl;

    cout << "Command line: ";
    for (int i = 0 ; i < argc ; i++){
        cout << argv[i] << " ";
    }
    cout << endl;
    cout << "Alice Assembler version " << version << "\nLast update: " << date << "\nAuthor: " << author << endl << endl;

    if(!clipp::parse(argc, argv, cli)) {
        if (!help){
            cout << "Could not parse the arguments" << endl;
            cout << clipp::make_man_page(cli, argv[0]);
            exit(1);
        }
        else{
            cout << "Help: " << endl;
            cout << clipp::make_man_page(cli, argv[0]);

            string command_bcalm = path_to_bcalm + " --help 2> trash.log > trash.log";
            auto bcalm_ok = system(command_bcalm.c_str());

            if (bcalm_ok == 0){
                exit(0);
            }
            else{
                cout << "Missing dependency: bcalm" << endl;
                exit(1);
            }
        }
    }

    if (order % 2 == 0){
        cerr << "WARNING: order (-l) must be odd, changing l to " << order-1 << "\n";
        order = order-1;
    }

    if (assembler != "custom" && assembler != "hifiasm" && assembler != "spades" && assembler != "gatb-minia" && assembler != "raven" && assembler != "megahit"){
        cerr << "ERROR: assembler must be bcalm or hifiasm or spades or gatb-mina or raven or flye \n";
        exit(1);
    }
    vector<int> kmer_sizes_vector;
    if (assembler == "custom"){
        std::stringstream ss(kmer_sizes);
        string item;
        while (std::getline(ss, item, ',')) {
            try {
            int k = std::stoi(item);
            if (!kmer_sizes_vector.empty() && k <= kmer_sizes_vector.back()) {
                cerr << "ERROR: kmer-sizes must be in increasing order\n";
                exit(1);
            }
            kmer_sizes_vector.push_back(k);
            } catch (const std::invalid_argument& e) {
            cerr << "ERROR: kmer-sizes must be a list of integers\n";
            exit(1);
            }
        }
        //if the last value <31, add 31 to do a proper expansion (>=km)
        if (kmer_sizes_vector.size()== 0 || kmer_sizes_vector[kmer_sizes_vector.size()-1] < 31){
            cerr << "WARNING: the kmer_sizes must go at least to 31, adding 31 at the end of your list" << endl;
            kmer_sizes_vector.push_back(31);
        }
    }

    //make sure the output folder ends with a /
    if (output_folder[output_folder.size()-1] != '/'){
        output_folder += "/";
    }
    string tmp_folder = output_folder + "tmp/";

    //record time now to measure the time of the whole process
    auto start = std::chrono::high_resolution_clock::now();

    //check if the output folder exists
    string command = "mkdir -p " + output_folder + " 2> trash.log";
    auto res = system(command.c_str());

    //create the tmp folder
    command = "mkdir -p " + tmp_folder + " 2> trash.log";
    res = system(command.c_str());
    
    string path_src = argv[0];
    path_src = path_src.substr(0, path_src.find_last_of("/")); //strip the /aliceasm
    path_src = path_src.substr(0, path_src.find_last_of("/")); //strip the /build

    std::string path_convertToGFA = "python3 " + path_src + "/bcalm/scripts/convertToGFA.py";
    string path_graphunzip = path_src + "/build/graphunzip";

    string path_total = argv[0];
    string path_fastg2gfa = path_total.substr(0, path_total.find_last_of("/"))+ "/fastg2gfa";

    check_dependencies(assembler, path_to_bcalm, path_to_hifiasm, path_to_spades, path_to_minia, path_to_raven, path_to_flye, path_to_minimap2, path_to_miniasm, path_to_minipolish, path_to_megahit, path_fastg2gfa, path_convertToGFA, path_graphunzip, path_src);

    // Check if the input file is gzipped
    if (input_file.substr(input_file.find_last_of('.') + 1) == "gz") {
        string unzipped_file = input_file.substr(input_file.find_last_of('/') + 1, input_file.find_last_of('.') - input_file.find_last_of('/') - 1);
        string command = "gunzip -c " + input_file + " > " + tmp_folder + unzipped_file;
        auto ok = system(command.c_str());
        if (ok != 0) {
            cerr << "ERROR: Failed to unzip the gzipped input file\n";
            exit(1);
        }
        input_file = tmp_folder + unzipped_file;
    }
    //if the input file is a fastq file, convert it to fasta
    if (input_file.substr(input_file.find_last_of('.')+1) == "fastq" || input_file.substr(input_file.find_last_of('.')+1) == "fq"){
        string fasta_file = input_file.substr(input_file.find_last_of('/') + 1, input_file.find_last_of('.') - input_file.find_last_of('/') - 1) + ".fasta";
        string command = "sed -n '1~4s/^@/>/p;2~4p' " + input_file + " > " + tmp_folder + fasta_file;
        auto ok = system(command.c_str());
        if (ok != 0){
            cerr << "ERROR: fastq_to_fasta failed\n";
            exit(1);
        }
        input_file = tmp_folder+fasta_file;
    }

    if (single_genome){
        contiguity = true;
    }

    string compressed_file = tmp_folder+"compressed.fa";
    string sampled_file = tmp_folder+"sampled.fa";
    string merged_gfa = tmp_folder+"bcalm.unitigs.shaved.merged.gfa";
    int km = 31; //size of the kmer used to do the expansion. Must be >21

    cout << "==== Step 1: MSR compression of the reads ====" << endl;
    
    auto time_start = std::chrono::high_resolution_clock::now();
    reduce(input_file, compressed_file, order, compression, num_threads, homopolymer_compression);
    auto time_reduced = std::chrono::high_resolution_clock::now();

    cout << "Done compressing reads, the compressed reads are in " << compressed_file << "\n" << endl;

    cout << "==== Step 2: Assembly of the compressed reads with " + assembler + " ====" << endl;
    string compressed_assembly = tmp_folder+"assembly_compressed.gfa";
    if (assembler == "custom"){
        assembly_custom(compressed_file, min_abundance, contiguity, (int) 20000/compression, tmp_folder, num_threads, compressed_assembly, kmer_sizes_vector, single_genome, path_to_bcalm, path_convertToGFA, path_graphunzip);
    }
    else if (assembler == "hifiasm"){
        assembly_hifiasm(compressed_file, tmp_folder, num_threads, compressed_assembly, path_to_hifiasm, assembler_parameters);
    }
    else if (assembler == "spades"){
        assembly_spades(compressed_file, tmp_folder, num_threads, compressed_assembly, path_to_spades, assembler_parameters);
    }
    else if (assembler == "gatb-minia"){
        assembly_minia(compressed_file, tmp_folder, num_threads, compressed_assembly, path_to_minia, path_convertToGFA, assembler_parameters);
    }
    else if (assembler == "raven"){
        assembly_raven(compressed_file, tmp_folder, num_threads, compressed_assembly, path_to_raven, assembler_parameters);
    }
    else if (assembler == "flye"){
        assembly_flye(compressed_file, tmp_folder, num_threads, compressed_assembly, path_to_flye, assembler_parameters);
    }
    else if (assembler == "miniasm"){
        assembly_miniasm(compressed_file, tmp_folder, num_threads, compressed_assembly, path_to_miniasm, path_to_minimap2, path_to_minipolish, assembler_parameters);
    }
    else if (assembler == "megahit"){
        assembly_megahit(compressed_file, tmp_folder, num_threads, compressed_assembly, path_to_megahit, path_fastg2gfa, assembler_parameters);
    }

    auto time_assembled = std::chrono::high_resolution_clock::now();

    cout << "==== Step 3: Inflating back the assembly to non-compressed space ====\n";


    //now let's parse the gfa file and decompress it
    time_t now2 = time(0);
    tm *ltm2 = localtime(&now2);
    cout << " - Listing the kmers needed for the expansion [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;

    unordered_map<uint64_t, pair<unsigned long long, unsigned long long>> kmers;
    std::unordered_set<uint64_t> central_kmers_needed;
    std::unordered_set<uint64_t> full_kmers_needed;
    string decompressed_assembly = tmp_folder+"assembly_decompressed.gfa";
    string central_kmer_file = tmp_folder+"central_kmers.txt";
    string full_kmer_file = tmp_folder+"full_kmers.txt";

    //list_kmers_needed_for_expansion(compressed_assembly, km, kmers_needed);
    expand_or_list_kmers_needed_for_expansion("index", compressed_assembly, km, central_kmers_needed, full_kmers_needed, central_kmer_file, full_kmer_file, kmers, decompressed_assembly);

    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << " - Parsing the reads to map compressed kmers with uncompressed sequences [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    go_through_the_reads_again_and_index_interesting_kmers(input_file, compressed_assembly, order, compression, km, central_kmers_needed, full_kmers_needed, kmers, central_kmer_file, full_kmer_file, num_threads, homopolymer_compression);


    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << " - Reconstructing the uncompressed assembly [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    //expand(compressed_assembly, decompressed_assembly, km, kmer_file, kmers);
    expand_or_list_kmers_needed_for_expansion("expand", compressed_assembly, km, central_kmers_needed, full_kmers_needed, central_kmer_file, full_kmer_file, kmers, decompressed_assembly);

    string output_file = output_folder + "assembly.gfa";
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << " - Computing the exact overlaps between the contigs [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    compute_exact_CIGARs(decompressed_assembly, output_file, 150*compression, 70*compression);

    //convert to fasta
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << " - Converting the assembly to fasta [" << 1+ ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    gfa_to_fasta(output_file, output_file.substr(0, output_file.find_last_of('.')) + ".fasta");

    //clean the tmp folder if the user wants
    if (clean){
        command = "rm -r " + tmp_folder + " 2> trash.log";
        auto res = system(command.c_str());
    }
    auto time_end = std::chrono::high_resolution_clock::now();

    cout << "\nDone, the final assembly is in " << output_file << "\n" << endl;
    cout << "Timing:\n";
    cout << "Compression: " << std::chrono::duration_cast<std::chrono::seconds>(time_reduced - time_start).count() << "s\n";
    cout << "Assembly: " << std::chrono::duration_cast<std::chrono::seconds>(time_assembled - time_reduced).count() << "s\n";
    cout << "Decompression: " << std::chrono::duration_cast<std::chrono::seconds>(time_end - time_assembled).count() << "s\n";
    cout << "Total time: " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count() << "s\n";

    if (test_ref_gfa != ""){
        cout << "Comparing the result with the reference " << test_ref_gfa << endl;
        test_assembly(output_file, test_ref_gfa);
    }

    // Remove the trash.log file if it exists
    std::remove("trash.log");

    return 0;
}

