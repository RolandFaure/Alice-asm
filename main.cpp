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

string version = "0.6.2";
string date = "2024-07-17";
string author = "Roland Faure";

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


    // string command_minia = path_minia + " --help 2> trash.log > trash.log";
    // auto minia_ok = system(command_minia.c_str());
    // cout << "trying command line " << command_minia << " " << minia_ok << " " << minia_ok2 << endl;

    int convertToGFA_ok = 0;
    int graphunzip_ok = 0;
    if (assembler == "bcalm"){
        path_convertToGFA = "python " + path_src + "/bcalm/scripts/convertToGFA.py";
        convertToGFA_ok = system((path_convertToGFA+" --help >trash.log 2>trash.log").c_str());
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

        auto graphunzip_ok = system((path_graphunzip + " --help >trash.log 2>trash.log ").c_str());
        if (graphunzip_ok != 0){
            path_graphunzip = "graphunzip";
            graphunzip_ok = system(path_graphunzip.c_str());
            if (graphunzip_ok != 0){
                cerr << "ERROR: graphunzip not found, problem in the installation, error code 322.\n";
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
    std::cout << "-------------------------------" << std::endl;


    if ((bcalm_ok != 0 && assembler == "bcalm") || 
        (hifiasm_ok != 0 && assembler == "hifiasm") || 
        (spades_ok != 0 && assembler == "spades") || 
        // (minia_ok != 0 && assembler == "gatb-minia") ||
        (raven_ok != 0 && assembler == "raven") ||
        (flye_ok != 0 && assembler == "flye") ||
        ((minimap_ok != 0 || miniasm_ok != 0 || minipolish_ok != 0) && assembler == "miniasm") ||
        (megahit_ok != 0 && assembler == "megahit" && fastg2gfa_ok != 0) ||
        (convertToGFA_ok != 0 || graphunzip_ok != 0)){
        std::cout << "Error: some dependencies are missing. Please install them or provide a valid path with the options." << std::endl;
        exit(1);
    }

}

int main(int argc, char** argv)
{

    string shaved_and_popped_gfa2 = "bcalm.unitigs.shaved.popped.gfa";
    string merged_gfa2 = "bcalm.unitigs.shaved.merged.gfa";
    int size_longest_read = 1000;
    pop_bubbles(merged_gfa2, size_longest_read, shaved_and_popped_gfa2);
    unordered_map<string, int> segments_IDs;
    vector<Segment> segments;
    vector<Segment> merged_segments;
    load_GFA(shaved_and_popped_gfa2, segments, segments_IDs);
    string shaved_and_popped_merged = "bcalm.unitigs.shaved.popped.merged.gfa";
    merge_adjacent_contigs(segments, merged_segments, shaved_and_popped_gfa2);
    output_graph(shaved_and_popped_merged, shaved_and_popped_gfa2, merged_segments);

    exit(1);

    //use clipp to parse the command line
    bool help = false;
    string input_file, output_folder;
    string assembler = "bcalm";
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
    bool rescue = false;
    bool contiguity = true;
    int min_abundance = 5;
    int order = 201;
    int compression = 20;
    int num_threads = 1;
    bool no_hpc = false;
    auto cli = (
        //input/output option
        clipp::required("-r", "--reads").doc("input file (fasta/q)") & clipp::opt_value("r",input_file),
        clipp::required("-o", "--output").doc("output folder") & clipp::opt_value("o",output_folder),

        //Performance options
        clipp::option("-t", "--threads").doc("number of threads [1]") & clipp::opt_value("t", num_threads),

        //Compression options
        clipp::option("-l", "--order").doc("order of MSR compression (odd) [201]") & clipp::opt_value("o", order),
        clipp::option("-c", "--compression").doc("compression factor [20]") & clipp::opt_value("c", compression),
        clipp::option("-H", "--no-hpc").set(no_hpc).doc("turn off homopolymer compression"),

        //Assembly options for the custom assembler
        clipp::option("-m", "--min-abundance").doc("minimum abundance of kmer to consider solid (for custom Alice assembler) [5]") & clipp::opt_value("m", min_abundance),
        clipp::option("--contiguity").set(contiguity).doc("Favor contiguity over recovery of rare strains [off]"),

        //Other assemblers options
        clipp::option("-a", "--assembler").doc("assembler to use {bcalm, hifiasm, spades, raven, gatb-minia, megahit} [bcalm]") & clipp::opt_value("a", assembler),
        clipp::option("--parameters").doc("extra parameters to pass to the assembler (between quotation marks) [\"\"]") & clipp::opt_value("p", assembler_parameters),
        clipp::option("--bcalm").doc("path to bcalm [bcalm]") & clipp::opt_value("b", path_to_bcalm),
        clipp::option("--hifiasm_meta").doc("path to hifiasm_meta [hifiasm_meta]") & clipp::opt_value("h", path_to_hifiasm),
        clipp::option("--spades").doc("path to spades [spades.py]") & clipp::opt_value("s", path_to_spades),
        clipp::option("--raven").doc("path to raven [raven]") & clipp::opt_value("r", path_to_bcalm),
        // clipp::option("--flye").doc("path to flye [flye]") & clipp::opt_value("f", path_to_flye), //flye does not work well with compressed reads
        clipp::option("--gatb-minia").doc("path to gatb-minia [gatb]") & clipp::opt_value("g", path_to_minia),
        clipp::option("--megahit").doc("path to megahit [megahit]") & clipp::opt_value("m", path_to_megahit),
        // clipp::option("--miniasm").doc("path to miniasm [miniasm]") & clipp::opt_value("m", path_to_miniasm), //we did not manage to make miniasm work
        // clipp::option("--minimap2").doc("path to minimap2 [minimap2]") & clipp::opt_value("m", path_to_minimap2),
        // clipp::option("--minipolish").doc("path to minipolish [minipolish]") & clipp::opt_value("r", path_to_minipolish),
        

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
        cout << "Could not parse the arguments" << endl;
        cout << clipp::make_man_page(cli, argv[0]);
        exit(1);
    }

    if (order % 2 == 0){
        cerr << "WARNING: order (-l) must be odd, changing l to " << order-1 << "\n";
        order = order-1;
    }
    int context_length = (order-1)/2;

    if (assembler != "bcalm" && assembler != "hifiasm" && assembler != "spades" && assembler != "gatb-minia" && assembler != "raven" && assembler != "megahit"){
        cerr << "ERROR: assembler must be bcalm or hifiasm or spades or gatb-mina or raven or flye \n";
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
    string command = "mkdir -p " + output_folder + " 2> trash.log";
    system(command.c_str());

    //create the tmp folder
    command = "mkdir -p " + tmp_folder + " 2> trash.log";
    system(command.c_str());
    

    string path_src = argv[0];
    path_src = path_src.substr(0, path_src.find_last_of("/")); //strip the /reduce
    path_src = path_src.substr(0, path_src.find_last_of("/")); //strip the /build

    std::string path_convertToGFA = "python " + path_src + "/bcalm/scripts/convertToGFA.py";
    string path_graphunzip = path_src + "/build/graphunzip";

    string path_total = argv[0];
    string path_fastg2gfa = path_total.substr(0, path_total.find_last_of("/"))+ "/fastg2gfa";

    check_dependencies(assembler, path_to_bcalm, path_to_hifiasm, path_to_spades, path_to_minia, path_to_raven, path_to_flye, path_to_minimap2, path_to_miniasm, path_to_minipolish, path_to_megahit, path_fastg2gfa, path_convertToGFA, path_graphunzip, path_src);

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

    cout << "==== Step 1: MSR compression of the reads ====" << endl;
    
    reduce(input_file, compressed_file, context_length, compression, num_threads, homopolymer_compression);

    cout << "Done compressing reads, the compressed reads are in " << compressed_file << "\n" << endl;

    cout << "==== Step 2: Assembly of the compressed reads with " + assembler + " ====" << endl;
    string compressed_assembly = tmp_folder+"assembly_compressed.gfa";
    if (assembler == "bcalm"){
        assembly_bcalm(compressed_file, min_abundance, contiguity, (int) 20000/compression, tmp_folder, num_threads, compressed_assembly, path_to_bcalm, path_convertToGFA, path_graphunzip, assembler_parameters);
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

    cout << "==== Step 3: Inflating back the assembly to non-compressed space ====\n";

    //now let's parse the gfa file and decompress it
    cout << " - Parsing the reads to map compressed kmers with uncompressed sequences\n";

    unordered_map<string, pair<unsigned long long, unsigned long long>> kmers;
    std::set<string> kmers_needed;
    string decompressed_assembly = tmp_folder+"assembly_decompressed.gfa";
    list_kmers_needed_for_expansion(compressed_assembly, km, kmers_needed);
    string kmer_file = tmp_folder+"kmers.txt";
    go_through_the_reads_again_and_index_interesting_kmers(input_file, compressed_assembly, context_length, compression, km, kmers_needed, kmers, kmer_file, num_threads, homopolymer_compression);
    cout << " - Reconstructing the uncompressed assembly" << endl;
    expand(compressed_assembly, decompressed_assembly, km, kmer_file, kmers);

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

