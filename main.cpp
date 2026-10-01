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
#include <algorithm>
#include <omp.h> //for efficient parallelization
#include <set>
#include <filesystem>

#include "reduce_and_expand.h"
#include "basic_graph_manipulation.h"
#include "robin_hood.h"
#include "assembly.h"
#include "clipp.h"
#include "bluntify.h"

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

string version = "0.8.0";
string date = "2026-10-01";
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

//returns true if the command runs and exits with the expected code (output discarded)
static bool runs(const string& command, int expected_code = 0){
    return system((command + " > /dev/null 2>&1").c_str()) == expected_code;
}

void check_dependencies(string assembler, string path_bcalm, string path_hifiasm, string path_spades, string path_minia, string path_raven, string path_megahit, string path_fastg2gfa,
    string &path_convertToGFA, string &path_graphunzip){

    bool python3_ok = runs("python3 --version");

    bool convertToGFA_ok = runs(path_convertToGFA + " -h");
    if (!convertToGFA_ok) {
        string bad_path = path_convertToGFA;
        string which_convertToGFA = exec("which convertToGFA.py");
        path_convertToGFA = "python3 " + shell_quote(which_convertToGFA);
        convertToGFA_ok = which_convertToGFA != "" && runs(path_convertToGFA + " -h");
        if (!convertToGFA_ok) {
            cerr << "ERROR: convertToGFA.py not found, problem in the installation, error code 321.\n";
            cout << "tried " << endl << bad_path << endl << path_convertToGFA << endl;
            exit(1);
        }
    }

    if (assembler == "custom" && !runs(path_graphunzip + " --help")){
        path_graphunzip = "graphunzip";
        if (!runs(path_graphunzip + " --help")){
            cerr << "ERROR: graphunzip not found, problem in the installation, error code 322.\n";
            exit(1);
        }
    }

    //only check the tools needed by the chosen assembler
    vector<pair<string, bool>> tools;
    if (assembler == "custom")
        tools = {{"bcalm", runs(path_bcalm + " --help")}};
    else if (assembler == "hifiasm")
        tools = {{"hifiasm", runs(path_hifiasm + " -h")}};
    else if (assembler == "spades")
        tools = {{"spades", runs(path_spades + " --help")}};
    else if (assembler == "gatb-minia")
        tools = {{"gatb-minia", runs(path_minia + " --help")}};
    else if (assembler == "raven")
        tools = {{"raven", runs(path_raven + " --help")}};
    else if (assembler == "megahit")
        tools = {{"megahit", runs(path_megahit + " --version")}, {"fastg2gfa", runs(path_fastg2gfa, 256)}};

    auto row = [](string name, bool ok){
        name.resize(15, ' ');
        std::cout << "|    " << name << "|   " << (ok ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    };
    std::cout << "_______________________________" << std::endl;
    std::cout << "|    Dependency     |  Found  |" << std::endl;
    std::cout << "|-------------------|---------|" << std::endl;
    row("python3", python3_ok);
    bool all_ok = python3_ok;
    for (auto& t : tools){
        row(t.first, t.second);
        all_ok = all_ok && t.second;
    }
    std::cout << "-------------------------------" << std::endl;

    if (!all_ok){
        std::cout << "Error: some dependencies are missing." << std::endl;
        exit(1);
    }
}

int main(int argc, char** argv)
{
    //use clipp to parse the command line
    bool help = false;
    string input_file, output_folder;
    string assembler = "custom";
    string path_to_bcalm = "bcalm";
    string path_to_hifiasm = "hifiasm_meta";
    string path_to_spades = "spades.py";
    string path_to_minia = "gatb";
    string path_to_raven = "raven";
    string path_to_megahit = "megahit";
    string assembler_parameters = "";
    bool contiguity = false;
    bool single_genome= false;
    int min_abundance = 5;
    string kmer_sizes = "21,31,61,101,191";
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
        clipp::option("-H", "--no-hpc").set(no_hpc).doc("turn off homopolymer compression"),

        //Assembly options for the custom assembler
        clipp::option("-m", "--min-abundance").doc("minimum abundance of kmer to consider solid - RECOMMENDED to set to coverage/2 if single-genome [5]") & clipp::opt_value("m", min_abundance),
        clipp::option("-k", "--kmer-sizes").doc("comma-separated increasing sizes of k for assembly, must go at least to 31 [21,31,61,101,191]") & clipp::opt_value("k", kmer_sizes),
        clipp::option("--single-genome").set(single_genome).doc("Switch on if assembling a single genome"),
        clipp::option("--contiguity").set(contiguity).doc("Favors contiguity by popping bubbles in the gfa graph [off]"),

        //Other assemblers options
        clipp::option("-a", "--assembler").doc("assembler to use {custom, spades} [custom]") & clipp::opt_value("a", assembler),
        
        // clipp::option("--parameters").doc("extra parameters to pass to the assembler (between quotation marks) [\"\"]") & clipp::opt_value("p", assembler_parameters),
        clipp::option("--bcalm").doc("path to bcalm [bcalm]") & clipp::opt_value("b", path_to_bcalm),
        // clipp::option("--hifiasm_meta").doc("path to hifiasm_meta [hifiasm_meta]") & clipp::opt_value("h", path_to_hifiasm),
        clipp::option("--spades").doc("path to spades [spades.py]") & clipp::opt_value("s", path_to_spades),
        // clipp::option("--raven").doc("path to raven [raven]") & clipp::opt_value("r", path_to_bcalm),
        // clipp::option("--gatb-minia").doc("path to gatb-minia [gatb]") & clipp::opt_value("g", path_to_minia),
        // clipp::option("--megahit").doc("path to megahit [megahit]") & clipp::opt_value("m", path_to_megahit),
        
        //Other options
        clipp::option("--clean").set(clean).doc("remove the tmp folder at the end [off]"),
        clipp::option("-v", "--version").call([]{ std::cout << "version " << version << "\nLast update: " << date << "\nAuthor: " << author << std::endl; exit(0); }).doc("print version and exit"),
        clipp::option("-h", "--help").set(help).doc("print this help message and exit")
    );


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

            if (runs(path_to_bcalm + " --help")){
                exit(0);
            }
            else{
                cout << "Missing dependency: bcalm" << endl;
                exit(1);
            }
        }
    }

    bool homopolymer_compression = !no_hpc; //must be read after parsing the command line

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
    if (output_folder.empty()){
        cerr << "ERROR: the output folder (-o) cannot be empty\n";
        exit(1);
    }
    if (output_folder[output_folder.size()-1] != '/'){
        output_folder += "/";
    }
    string tmp_folder = output_folder + "tmp/";

    //record time now to measure the time of the whole process
    auto start = std::chrono::high_resolution_clock::now();

    //create the output and tmp folders
    std::error_code error_creating_folder;
    std::filesystem::create_directories(tmp_folder, error_creating_folder);
    if (error_creating_folder){
        cerr << "ERROR: could not create the output folder " << tmp_folder << ": " << error_creating_folder.message() << "\n";
        exit(1);
    }
    
    string path_src = argv[0];
    path_src = path_src.substr(0, path_src.find_last_of("/")); //strip the /aliceasm
    path_src = path_src.substr(0, path_src.find_last_of("/")); //strip the /build

    std::string path_convertToGFA = "python3 " + shell_quote(path_src + "/bcalm/scripts/convertToGFA.py");
    string path_graphunzip = shell_quote(path_src + "/build/graphunzip");

    string path_total = argv[0];
    string path_fastg2gfa = shell_quote(path_total.substr(0, path_total.find_last_of("/"))+ "/fastg2gfa");

    check_dependencies(assembler, path_to_bcalm, path_to_hifiasm, path_to_spades, path_to_minia, path_to_raven, path_to_megahit, path_fastg2gfa, path_convertToGFA, path_graphunzip);

    //gzipped, FASTQ, multi-line or lowercase inputs are converted to single-line uppercase FASTA
    input_file = prepare_reads(input_file, tmp_folder);

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
    int assembly_k = kmer_sizes_vector.empty() ? 191 : kmer_sizes_vector.back(); //k of the final compressed graph, bounds the overlaps between contigs
    if (assembler == "custom"){
        assembly_k = assembly_custom(compressed_file, min_abundance, tmp_folder, num_threads, compressed_assembly, kmer_sizes_vector, single_genome, path_to_bcalm, path_convertToGFA, path_graphunzip, contiguity);
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
    else if (assembler == "megahit"){
        assembly_megahit(compressed_file, tmp_folder, num_threads, compressed_assembly, path_to_megahit, path_fastg2gfa, assembler_parameters);
    }

    auto time_assembled = std::chrono::high_resolution_clock::now();

    cout << "==== Step 3: Inflating back the assembly to non-compressed space ====\n";


    //now let's parse the gfa file and decompress it
    time_t now2 = time(0);
    tm *ltm2 = localtime(&now2);
    cout << " - Listing the kmers needed for the expansion [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;

    unordered_map<uint64_t, pair<unsigned long long, unsigned long long>> kmers;
    std::vector<uint64_t> central_kmers_needed;
    std::vector<uint64_t> full_kmers_needed;
    string decompressed_assembly = tmp_folder+"assembly_decompressed.gfa";
    string central_kmer_file = tmp_folder+"central_kmers.txt";
    string full_kmer_file = tmp_folder+"full_kmers.txt";

    //list_kmers_needed_for_expansion(compressed_assembly, km, kmers_needed);
    expand_or_list_kmers_needed_for_expansion("index", compressed_assembly, km, compression, central_kmers_needed, full_kmers_needed, central_kmer_file, full_kmer_file, kmers, decompressed_assembly);

    // Sort and deduplicate the kmers needed
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << " - Sorting kmer vectors for efficient lookup [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    std::sort(central_kmers_needed.begin(), central_kmers_needed.end());
    std::sort(full_kmers_needed.begin(), full_kmers_needed.end());
    
    // Remove duplicates to optimize memory and search performance
    central_kmers_needed.erase(std::unique(central_kmers_needed.begin(), central_kmers_needed.end()), central_kmers_needed.end());
    full_kmers_needed.erase(std::unique(full_kmers_needed.begin(), full_kmers_needed.end()), full_kmers_needed.end());

    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << " - Parsing the reads to map compressed kmers with uncompressed sequences [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    go_through_the_reads_again_and_index_interesting_kmers(input_file, compressed_assembly, order, compression, km, central_kmers_needed, full_kmers_needed, kmers, central_kmer_file, full_kmer_file, num_threads, homopolymer_compression);


    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << " - Reconstructing the uncompressed assembly [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    //expand(compressed_assembly, decompressed_assembly, km, kmer_file, kmers);
    expand_or_list_kmers_needed_for_expansion("expand", compressed_assembly, km, compression, central_kmers_needed, full_kmers_needed, central_kmer_file, full_kmer_file, kmers, decompressed_assembly);

    string output_file = output_folder + "assembly.gfa";
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << " - Computing the exact overlaps between the contigs [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    long long compressed_length = total_sequence_length(compressed_assembly);
    double bases_per_compressed_base = compressed_length > 0 ? total_sequence_length(decompressed_assembly) / (double) compressed_length : compression;
    compute_exact_CIGARs(decompressed_assembly, output_file, assembly_k * 2 * compression, assembly_k * 1 * compression, bases_per_compressed_base, num_threads);

    if (single_genome){
        // bluntify the graph for single-genome assemblies
        bluntify(output_file, output_file, assembly_k * compression * 0.8, tmp_folder);
    }

    //convert to fasta
    now2 = time(0);
    ltm2 = localtime(&now2);
    cout << " - Converting the assembly to fasta [" << ltm2->tm_mday << "/" << 1 + ltm2->tm_mon << "/" << 1900 + ltm2->tm_year << " " << ltm2->tm_hour << ":" << ltm2->tm_min << ":" << ltm2->tm_sec << "]" << endl;
    gfa_to_fasta(output_file, output_file.substr(0, output_file.find_last_of('.')) + ".fasta");


    //clean the tmp folder if the user wants
    if (clean){
        std::filesystem::remove_all(tmp_folder);
    }
    auto time_end = std::chrono::high_resolution_clock::now();

    cout << "\nDone, the final assembly is in " << output_file << "\n" << endl;
    cout << "Timing:\n";
    cout << "Compression: " << std::chrono::duration_cast<std::chrono::seconds>(time_reduced - time_start).count() << "s\n";
    cout << "Assembly: " << std::chrono::duration_cast<std::chrono::seconds>(time_assembled - time_reduced).count() << "s\n";
    cout << "Decompression: " << std::chrono::duration_cast<std::chrono::seconds>(time_end - time_assembled).count() << "s\n";
    cout << "Total time: " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count() << "s\n";

    return 0;
}

