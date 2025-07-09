# Alice

Alice is a genomic and metagenomic assembler designed for HiFi read (<0.1% error rate). 


It is named after the character of Lewis Caroll's _Alice in Wonderland_, and more precisely to the "drink-me potion" and "eat-me cake", which make Alice respectively shrink and grow. The idea of the Alice assembler is to shrink the input reads, perform assembly on shrunken data and inflate the obtained assembly back to normal size. The compression is a Mapping-friendly Sequence Reduction (MSR) of high order, which allows all base to be taken into account during compression. However, it increases error rate of the reads. It is thus recommended to run it on long and accurate reads.
![alice_compression](https://github.com/rolandfaure/Alice-asm/blob/master/alice_compression.png)

## Installation

The recommended way to install Alice is through conda
```
conda install -c bioconda aliceasm
```

To download and compile the Alice assembler from source
```
git clone https://github.com/RolandFaure/Alice-asm.git
cd Alice-asm
mkdir build && cd build
cmake ..
make
```

## Usage

```
/build/aliceasm --help
   _______ _                     _ _                                            _     _                               
  |__   __| |              /\   | (_)              /\                          | |   | |                              
     | |  | |__   ___     /  \  | |_  ___ ___     /  \   ___ ___  ___ _ __ ___ | |__ | | ___ _ __      ::     _.mnm._ 
     | |  | '_ \ / _ \   / /\ \ | | |/ __/ _ \   / /\ \ / __/ __|/ _ \ '_ ` _ \| '_ \| |/ _ \ '__|    :  :   ( _____ )
     | |  | | | |  __/  / ____ \| | | (_|  __/  / ____ \\__ \__ \  __/ | | | | | |_) | |  __/ |       :  :    |     | 
     |_|  |_| |_|\___| /_/    \_\_|_|\___\___| /_/    \_\___/___/\___|_| |_| |_|_.__/|_|\___|_|       :__:     `___/  

Command line: ./build/aliceasm --help 
Alice Assembler version 0.6.41
Last update: 2024-04-09
Author: Roland Faure

Help: 
SYNOPSIS
        ./build/aliceasm -r [<r>] -o [<o>] [-t [<t>]] [-l [<o>]] [-c [<c>]] [-H] [-m [<m>]] [-k
                         [<k>]] [--single-genome] [--bcalm [<b>]] [--clean] [--test [<t>]] [-v] [-h]

OPTIONS
        -r, --reads input file (fasta/q)
        -o, --output
                    output folder

        -t, --threads
                    number of threads [1]

        -l, --order order of MSR compression (odd) [101]
        -c, --compression
                    compression factor [20]

        -H, --no-hpc
                    turn off homopolymer compression

        -m, --min-abundance
                    minimum abundance of kmer to consider solid - RECOMMENDED to set to coverage/2
                    if single-genome [5]

        -k, --kmer-sizes
                    comma-separated increasing sizes of k for assembly, must go at least to 31
                    [17,31]

        --single-genome
                    Switch on if assembling a single genome

        --bcalm     path to bcalm [bcalm]
        --clean     remove the tmp folder at the end [off]
        -v, --version
                    print version and exit

        -h, --help  print this help message and exit
```

## Can I use Alice for both genomic and metagenomic use cases ?

Yes! If you are trying to assemble a single genome, it is recommended to use flag --single-genome and set -m to expected-coverage/2 in order to simplify the output graph.
