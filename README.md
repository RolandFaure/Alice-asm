# The Alice Assembler

The Alice assembler is named after the character of Lewis Caroll's _Alice in Wonderland_, and more precisely to the "drink-me potion" and "eat-me cake", which make Alice respectively shrink and grow. The idea of the Alice assembler is to shrink the input reads, perform plain De Bruijn Graph assembly on shrunken data and inflate the obtained assembly back to normal size. The compression is a Mapping-friendly Sequence Reduction (MSR) of high order, which allows all base to be taken into account during compression.
![alice_compression](https://github.com/rolandfaure/Alice-asm/blob/master/alice_compression.png)

## Installation

To have all the dependencies, it is advised to create first a conda environment with all the dependencies (bcalm, minimap2, openmp, g++, scipy, numpy)
```
conda create -n hairsplitter -c bioconda -c conda-forge bcalm openmp libgomp gxx gcc scipy numpy
conda activate hairsplitter
```

Then download and compile the Alice assembler
```
git clone https://github.com/RolandFaure/Alice-asm.git
cd Alice-asm
mkdir build && cd build
cmake ..
make
```

## Usage

```
SYNOPSIS
        aliceasm -r [<i>] -o [<o>] [-t [<t>]] [-m [<m>]] [-l [<o>]] [-c [<c>]]
            [--bcalm [<b>]] [--minimap2 [<m>]] [-v]

OPTIONS
        -r, --reads input file (fasta/q)
        -o, --output
                    output file (gfa)

        -t, --threads
                    number of threads [1]

        -m, --min-abundance
                    minimum abundance of kmer to consider solid [10]

        -l, --order order of MSR compression (odd) [201]
        -c, --compression
                    compression factor [20]

        --bcalm     path to bcalm [bcalm]
        --minimap2  path to minimap2 [minimap2]
        -v, --version
                    print version and exit

```
