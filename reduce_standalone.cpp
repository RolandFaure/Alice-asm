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

int main(int argc, char** argv)
{

    //use clipp to parse the command line
    bool help = false;
    string input_file, output_file;
    string assembler = "custom";
    bool rescue = false;
    bool contiguity = false;
    int min_abundance = 5;
    int order = 101;
    int compression = 20;
    int num_threads = 1;
    bool no_hpc = false;
    bool clean = false;
    auto cli = (
        //input/output option
        clipp::required("-s", "--seq").doc("input file (fasta/q)") & clipp::opt_value("r",input_file),
        clipp::required("-o", "--output").doc("output file") & clipp::opt_value("o",output_file),

        //Performance options
        clipp::option("-t", "--threads").doc("number of threads [1]") & clipp::opt_value("t", num_threads),

        //Compression options
        clipp::option("-l", "--order").doc("order of MSR compression (odd) [101]") & clipp::opt_value("o", order),
        clipp::option("-c", "--compression").doc("compression factor [20]") & clipp::opt_value("c", compression),
        clipp::option("-H", "--no-hpc").set(no_hpc).doc("turn off homopolymer compression")


    );

    bool homopolymer_compression = !no_hpc;

    if(!clipp::parse(argc, argv, cli)) {
        cout << "Could not parse the arguments" << endl;
        cout << clipp::make_man_page(cli, argv[0]);
        exit(1);
    }

    if (order % 2 == 0){
        cerr << "WARNING: order (-l) must be odd because of lousy software engineering, changing l to " << order-1 << "\n";
        order = order-1;
    }

    
    string compressed_file = output_file;

    int km = 31; //size of the kmer used to do the expansion. Must be >21

    auto time_start = std::chrono::high_resolution_clock::now();
    reduce(input_file, compressed_file, order, compression, num_threads, homopolymer_compression);
    auto time_reduced = std::chrono::high_resolution_clock::now();

    return 0;
}

