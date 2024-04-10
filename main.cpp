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

#include "reduce_and_expand.h"
#include "basic_graph_manipulation.h"
#include "robin_hood.h"
#include "assembly.h"
#include "clipp.h"

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

// ANSI escape codes for text color
#define RED_TEXT "\033[1;31m"
#define GREEN_TEXT "\033[1;32m"
#define RESET_TEXT "\033[0m"

string version = "0.2.0";
string date = "2024-04-10";
string author = "Roland Faure";

void check_dependencies(string assembler, string path_bcalm, string path_hifiasm, string path_spades, string path_minia, string path_raven, string &path_convertToGFA, string &path_graphunzip, string path_src){

    string command_bcalm = path_bcalm + " --help 2> /dev/null > /dev/null";
    auto bcalm_ok = system(command_bcalm.c_str());

    string command_hifiasm = path_hifiasm + " --help 2> /dev/null > /dev/null";
    auto hifiasm_ok = system(command_hifiasm.c_str());

    string command_spades = path_spades + " --help 2> /dev/null > /dev/null";
    auto spades_ok = system(command_spades.c_str());

    string command_raven = path_raven + " --help 2> /dev/null > /dev/null";
    auto raven_ok = system(command_raven.c_str());

    // string command_minia = path_minia + " --help 2> /dev/null > /dev/null";
    // auto minia_ok = system(command_minia.c_str());
    // cout << "trying command line " << command_minia << " " << minia_ok << " " << minia_ok2 << endl;

    int convertToGFA_ok = 0;
    int graphunzip_ok = 0;
    if (assembler == "bcalm"){
        path_convertToGFA = "python " + path_src + "/bcalm/scripts/convertToGFA.py";
        convertToGFA_ok = system((path_convertToGFA+" --help >/dev/null 2>/dev/null").c_str());
        if (convertToGFA_ok != 0){
            path_convertToGFA = "convertToGFA.py";
            convertToGFA_ok = system(path_convertToGFA.c_str());
            if (convertToGFA_ok != 0){
                cerr << "ERROR: convertToGFA.py not found, problem in the installation, error code 321. I was looking for convertToGFA.py at " << path_src << "/bcalm/scripts/convertToGFA.py but did not find it\n";
                cerr << "Command 1: " << "python " << path_src << "/bcalm/scripts/convertToGFA.py\n";
                cerr << "Command 2: " << "convertToGFA.py\n";
                exit(1);
            }
        }

        path_graphunzip = "python " + path_src + "/GraphUnzip/graphunzip.py";
        auto graphunzip_ok = system((path_graphunzip + " --help >/dev/null 2>/dev/null ").c_str());
        if (graphunzip_ok != 0){
            path_graphunzip = "graphunzip.py";
            graphunzip_ok = system(path_graphunzip.c_str());
            if (graphunzip_ok != 0){
                cerr << "ERROR: graphunzip.py not found, problem in the installation, error code 322. I was looking for graphunzip.py at " << path_src << "/GraphUnzip/graphunzip.py but did not find it\n";
                exit(1);
            }
        }
    }

    std::cout << "_______________________________" << std::endl;
    std::cout << "|    Dependency     |  Found  |" << std::endl;
    std::cout << "|-------------------|---------|" << std::endl;
    if (assembler == "bcalm")
        std::cout << "|    bcalm          |   " << (bcalm_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    else if (assembler == "hifiasm")
        std::cout << "|    hifiasm        |   " << (hifiasm_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    else if (assembler == "spades")
        std::cout << "|    spades         |   " << (spades_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    // else if (assembler == "gatb-minia")
    //     std::cout << "|    gatb-minia     |   " << (minia_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    else if (assembler == "raven")
        std::cout << "|    raven          |   " << (raven_ok == 0 ? GREEN_TEXT "Yes" : RED_TEXT "No ") << RESET_TEXT "   |" << std::endl;
    std::cout << "-------------------------------" << std::endl;

    if ((bcalm_ok != 0 && assembler == "bcalm") || 
        (hifiasm_ok != 0 && assembler == "hifiasm") || 
        (spades_ok != 0 && assembler == "spades") || 
        // (minia_ok != 0 && assembler == "gatb-minia") ||
        (raven_ok != 0 && assembler == "raven") ||
        convertToGFA_ok != 0 || graphunzip_ok != 0){
        std::cout << "Error: some dependencies are missing. Please install them or provide a valid path with the options." << std::endl;
        exit(1);
    }

}

int main(int argc, char** argv)
{

    //use clipp to parse the command line
    bool help = false;
    string input_file, output_folder;
    string assembler = "bcalm";
    string path_to_bcalm = "bcalm";
    string path_to_hifiasm = "hifiasm";
    string path_to_spades = "spades.py";
    string path_to_minia = "gatb";
    string path_to_raven = "raven";
    bool rescue = false;
    int min_abundance = 5;
    int order = 201;
    int compression = 20;
    int num_threads = 1;
    bool no_hpc = false;
    auto cli = (
        clipp::required("-r", "--reads").doc("input file (fasta/q)") & clipp::opt_value("r",input_file),
        clipp::required("-o", "--output").doc("output folder") & clipp::opt_value("o",output_folder),
        clipp::option("-t", "--threads").doc("number of threads [1]") & clipp::opt_value("t", num_threads),
        clipp::option("-a", "--assembler").doc("assembler to use {bcalm, hifiasm, spades, raven, gatb-minia} [bcalm]") & clipp::opt_value("a", assembler),
        clipp::option("-l", "--order").doc("order of MSR compression (odd) [201]") & clipp::opt_value("o", order),
        clipp::option("-c", "--compression").doc("compression factor [20]") & clipp::opt_value("c", compression),
        clipp::option("-H", "--no-hpc").set(no_hpc).doc("turn off homopolymer compression"),
        clipp::option("--bcalm").doc("path to bcalm [bcalm]") & clipp::opt_value("b", path_to_bcalm),
        clipp::option("--hifiasm").doc("paht to hifiasm [hifiasm]") & clipp::opt_value("h", path_to_hifiasm),
        clipp::option("--spades").doc("path to spades [spades.py]") & clipp::opt_value("s", path_to_spades),
        clipp::option("--raven").doc("path to raven [raven]") & clipp::opt_value("r", path_to_bcalm),
        clipp::option("--gatb-minia").doc("path to gatb-minia [gatb]") & clipp::opt_value("g", path_to_minia),
        clipp::option("-m", "--min-abundance").doc("minimum abundance of kmer to consider solid [5] (for kmer-based assemblers)") & clipp::opt_value("m", min_abundance),

        clipp::option("-v", "--version").call([]{ std::cout << "version " << version << "\nLast update: " << date << "\nAuthor: " << author << std::endl; exit(0); }).doc("print version and exit")

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


    //Show the bottle left, then the cake right
    cout << "          _______ _                     _ _                                            _     _              "           << "           " << endl;
    cout << "         |__   __| |              /\\   | (_)              /\\                          | |   | |             "         << "           " << endl;
    cout << " ::         | |  | |__   ___     /  \\  | |_  ___ ___     /  \\   ___ ___  ___ _ __ ___ | |__ | | ___ _ __    "         << "   _.mnm._ " << endl;
    cout << ":  :        | |  | '_ \\ / _ \\   / /\\ \\ | | |/ __/ _ \\   / /\\ \\ / __/ __|/ _ \\ '_ ` _ \\| '_ \\| |/ _ \\ '__|   "<< "  ( _____ )" << endl;
    cout << ":  :        | |  | | | |  __/  / ____ \\| | | (_|  __/  / ____ \\\\__ \\__ \\  __/ | | | | | |_) | |  __/ |      "      << "   |     | " << endl;
    cout << ":__:        |_|  |_| |_|\\___| /_/    \\_\\_|_|\\___\\___| /_/    \\_\\___/___/\\___|_| |_| |_|_.__/|_|\\___|_|      "  << "    `___/  " << endl;
    cout << endl;

    cout << "Command line: ";
    for (int i = 0 ; i < argc ; i++){
        cout << argv[i] << " ";
    }
    cout << endl;
    cout << "Alice Assembler version " << version << "\nLast update: " << date << "\nAuthor: " << author << endl << endl;

    if(!clipp::parse(argc, argv, cli)) {
        cout << "Could not parse the arguments" << endl;
        cout << clipp::make_man_page(cli, argv[0]);
        exit(1);
    }

    if (order % 2 == 0){
        cerr << "ERROR: order must be odd\n";
        exit(1);
    }
    int context_length = (order-1)/2;

    if (assembler != "bcalm" && assembler != "hifiasm" && assembler != "spades" && assembler != "gatb-minia" && assembler != "raven"){
        cerr << "ERROR: assembler must be bcalm or hifiasm or spades\n";
        exit(1);
    }

    //make sure the output folder ends with a /
    if (output_folder[output_folder.size()-1] != '/'){
        output_folder += "/";
    }
    string tmp_folder = output_folder + "tmp/";

    //record time now to measure the time of the whole process
    auto start = std::chrono::high_resolution_clock::now();

    //check if the output folder exists
    string command = "mkdir -p " + output_folder + " 2> /dev/null";
    system(command.c_str());

    //create the tmp folder
    command = "mkdir -p " + tmp_folder + " 2> /dev/null";
    system(command.c_str());
    

    string path_src = argv[0];
    path_src = path_src.substr(0, path_src.find_last_of("/")); //strip the /reduce
    path_src = path_src.substr(0, path_src.find_last_of("/")); //strip the /build

    std::string path_convertToGFA = "python " + path_src + "/bcalm/scripts/convertToGFA.py";
    string path_graphunzip = "python " + path_src + "/GraphUnzip/graphunzip.py";

    check_dependencies(assembler, path_to_bcalm, path_to_hifiasm, path_to_spades, path_to_minia, path_to_raven, path_convertToGFA, path_graphunzip, path_src);

    //if the input file is a fastq file, convert it to fasta
    if (input_file.substr(input_file.find_last_of('.')+1) == "fastq" || input_file.substr(input_file.find_last_of('.')+1) == "fq"){
        string fasta_file = input_file.substr(0, input_file.find_last_of('.')) + ".fasta";
        string command = "sed -n '1~4s/^@/>/p;2~4p' " + input_file + " > " + tmp_folder + fasta_file;
        auto ok = system(command.c_str());
        if (ok != 0){
            cerr << "ERROR: fastq_to_fasta failed\n";
            exit(1);
        }
        input_file = tmp_folder+fasta_file;
    }

    string compressed_file = tmp_folder+"compressed.fa";
    string merged_gfa = tmp_folder+"bcalm.unitigs.shaved.merged.gfa";
    int km = 31; //size of the kmer used to do the expansion. Must be >21

    cout << "==== Step 1: MSR compression of the reads ====\n";
    
    reduce(input_file, compressed_file, context_length, compression, num_threads, homopolymer_compression);

    cout << "Done compressing reads, the compressed reads are in " << compressed_file << "\n" << endl;

    cout << "==== Step 2: Assembly of the compressed reads with " + assembler + " ====\n";
    string compressed_assembly = tmp_folder+"assembly_compressed.gfa";
    if (assembler == "bcalm"){
        assembly_bcalm(compressed_file, min_abundance, tmp_folder, num_threads, compressed_assembly, path_to_bcalm, path_convertToGFA, path_graphunzip);
    }
    else if (assembler == "hifiasm"){
        assembly_hifiasm(compressed_file, tmp_folder, num_threads, compressed_assembly, path_to_hifiasm);
    }
    else if (assembler == "spades"){
        assembly_spades(compressed_file, tmp_folder, num_threads, compressed_assembly, path_to_spades);
    }
    else if (assembler == "gatb-minia"){
        assembly_minia(compressed_file, tmp_folder, num_threads, compressed_assembly, path_to_minia, path_convertToGFA);
    }
    else if (assembler == "raven"){
        assembly_raven(compressed_file, tmp_folder, num_threads, compressed_assembly, path_to_raven);
    }

    cout << "==== Step 3: Inflating back the assembly to non-compressed space ====\n";

    //now let's parse the gfa file and decompress it
    cout << " - Parsing the reads to map compressed kmers with uncompressed sequences\n";
    //time the next function
    unordered_map<string, pair<string,string>> kmers;
    go_through_the_reads_again(input_file, compressed_assembly, context_length, compression, km, kmers, num_threads, homopolymer_compression);
    string decompressed_assembly = tmp_folder+"assembly_decompressed.gfa";
    cout << " - Reconstructing the uncompressed assembly" << endl;
    expand(compressed_assembly, decompressed_assembly, km, 70, kmers);

    //test compute_eact_cigar
    string output_file = output_folder + "assembly.gfa";
    cout << " - Computing the exact overlaps between the contigs\n";
    compute_exact_CIGARs(decompressed_assembly, output_file, 150*compression, 70*compression);

    //convert to fasta
    
    cout << " - Converting the assembly to fasta\n";
    gfa_to_fasta(output_file, output_file.substr(0, output_file.find_last_of('.')) + ".fasta");

    cout << "\nDone, the final assembly is in " << output_file << "\n" << endl;
    cout << "Total time: " << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start).count() << "s\n";

    return 0;
}

