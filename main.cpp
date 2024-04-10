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

string version = "0.2.0";
string date = "2024-04-10";
string author = "Roland Faure";


int main(int argc, char** argv)
{

    //use clipp to parse the command line
    bool help = false;
    string input_file, output_folder;
    string assembler = "bcalm";
    string path_to_bcalm = "bcalm";
    string path_to_minimap = "minimap2";
    bool rescue = false;
    int min_abundance = 5;
    int order = 201;
    int compression = 20;
    int num_threads = 1;
    auto cli = (
        clipp::required("-r", "--reads").doc("input file (fasta/q)") & clipp::opt_value("r",input_file),
        clipp::required("-o", "--output").doc("output folder") & clipp::opt_value("o",output_folder),
        clipp::option("-t", "--threads").doc("number of threads [1]") & clipp::opt_value("t", num_threads),
        clipp::option("-m", "--min-abundance").doc("minimum abundance of kmer to consider solid [5]") & clipp::opt_value("m", min_abundance),
        clipp::option("-l", "--order").doc("order of MSR compression (odd) [201]") & clipp::opt_value("o", order),
        clipp::option("-c", "--compression").doc("compression factor [20]") & clipp::opt_value("c", compression),
        clipp::option("--bcalm").doc("path to bcalm [bcalm]") & clipp::opt_value("b", path_to_bcalm),
        clipp::option("--minimap2").doc("path to minimap2 [minimap2]") & clipp::opt_value("m", path_to_minimap),
        clipp::option("-v", "--version").call([]{ std::cout << "version " << version << "\nLast update: " << date << "\nAuthor: " << author << std::endl; exit(0); }).doc("print version and exit")

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

    //The Alice Assembler ascii art


    //Show the bottle left, then the cake right
    cout << "          _______ _                     _ _                                            _     _              " << "" << endl;
    cout << "         |__   __| |              /\\   | (_)              /\\                          | |   | |             " << "" << endl;
    cout << " ::         | |  | |__   ___     /  \\  | |_  ___ ___     /  \\   ___ ___  ___ _ __ ___ | |__ | | ___ _ __    " << "   _.mnm._ " << endl;
    cout << ":  :        | |  | '_ \\ / _ \\   / /\\ \\ | | |/ __/ _ \\   / /\\ \\ / __/ __|/ _ \\ '_ ` _ \\| '_ \\| |/ _ \\ '__|   "<< "  ( _____ )" << endl;
    cout << ":  :        | |  | | | |  __/  / ____ \\| | | (_|  __/  / ____ \\\\__ \\__ \\  __/ | | | | | |_) | |  __/ |      " << "   |     |" << endl;
    cout << ":__:        |_|  |_| |_|\\___| /_/    \\_\\_|_|\\___\\___| /_/    \\_\\___/___/\\___|_| |_| |_|_.__/|_|\\___|_|      " << "    `___/  "<< endl;
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

    if (assembler != "bcalm" && assembler != "hifiasm"){
        cerr << "ERROR: assembler must be bcalm or hifiasm\n";
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

    //check if the path to bcalm is correct
    string command_bcalm = path_to_bcalm + " --help 2> /dev/null > /dev/null";
    auto bcalm_ok = system(command_bcalm.c_str());
    if (bcalm_ok != 0){
        cerr << "ERROR: bcalm not found using command line " << command_bcalm << "\n";
        cerr << "Please specify a valid path using --bcalm\n" << endl;
        exit(1);
    }
    //check if the path to minimap is correct
    command = path_to_minimap + " --help 2> /dev/null > /dev/null";
    auto ok = system(command.c_str());
    if (ok != 0){
        cerr << "ERROR: minimap2 not found using command line " << command << "\n";
        cerr << "Please specify a valid path using --minimap2\n" << endl;
        exit(1);
    }

    std::string path_convertToGFA = "python " + path_src + "/bcalm/scripts/convertToGFA.py";

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
    
    reduce(input_file, compressed_file, context_length, compression, num_threads);

    cout << "Done compressing reads, the compressed reads are in " << compressed_file << "\n" << endl;

    cout << "==== Step 2: Assembly of the compressed reads with ====\n";
    string compressed_assembly = tmp_folder+"assembly_compressed.gfa";
    if (assembler == "bcalm"){
        assembly_bcalm(compressed_file, min_abundance, tmp_folder, num_threads, compressed_assembly, path_to_bcalm, path_convertToGFA, path_src);
    }
    else if (assembler == "hifiasm"){
        assembly_hifiasm(compressed_file, tmp_folder, num_threads, compressed_assembly);
    }

    cout << "==== Step 3: Inflating back the assembly to non-compressed space ====\n";

    //now let's parse the gfa file and decompress it
    cout << " - Parsing the reads to map compressed kmers with uncompressed sequences\n";
    //time the next function
    unordered_map<string, pair<string,string>> kmers;
    go_through_the_reads_again(input_file, compressed_assembly, context_length, compression, km, kmers, num_threads);
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

