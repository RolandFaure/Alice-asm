#include <gatb/bcalm2/bcalm_algo.cpp>
#include <gatb/bcalm2/bglue_algo.cpp>
#include <gatb/bcalm2/ograph.cpp>
#include <gatb/debruijn/impl/LinkTigs.cpp>

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

/********************************************************************************/
namespace gatb { namespace core { namespace debruijn { namespace impl  {
/********************************************************************************/


template void bcalm2<64>(Storage* storage, 
        std::string prefix,
        int kmerSize, 
        int abundance, 
        int minSize, 
        int nb_threads, 
        int minimizer_type, 
        bool verbose
        );
template void bglue<64>(Storage* storage, 
        std::string prefix,
        int kmerSize, 
        int nb_glue_partitions, 
        int nb_threads, 
        bool all_abundance_counts,
        bool verbose
        );

template class graph3<64>; // graph3<span> switch  

template void link_tigs<64>
    (std::string unitigs_filename, int kmerSize, int nb_threads, uint64_t &nb_unitigs, bool verbose, bool edge_km_representation, bool renumber_unitigs = false);

template void link_unitigs_pass<64>(const std::string unitigs_filename, bool verbose, const int pass, const int kmerSize, bool edge_km_representation, const bool renumber_unitigs);


/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
