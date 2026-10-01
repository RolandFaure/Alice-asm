#ifndef GRAPHUNZIP_H
#define GRAPHUNZIP_H

#include <iostream>
#include <vector>
#include <string>
#include <mutex>
#include <set>
#include <string>
#include <fstream>
#include <sstream>

#include "robin_hood.h"
#include "basic_graph_manipulation.h"

// The Segment class and load_GFA / output_graph / merge_adjacent_contigs are declared in basic_graph_manipulation.h.
// graphunzip.cpp defines the graphunzip executable. It keeps all the read-path information (memory-mapped GAF paths,
// haploidy, consensuses, strong neighbors) in its own structures and does not use the Segment path members
// (add_neighbor, compute_consensuses, get_strong_neighbors_*, neighbors_*, consensus_*, is_haploid) anymore.

#endif
