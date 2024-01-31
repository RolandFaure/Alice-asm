#include <gatb/debruijn/impl/Simplifications.cpp>
#include <gatb/debruijn/impl/UnitigsConstructionAlgorithm.cpp>
#include <gatb/debruijn/impl/GraphUnitigs.cpp>

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

/********************************************************************************/
namespace gatb { namespace core { namespace debruijn { namespace impl  {
/********************************************************************************/

typedef boost::variant<GraphData<64> > GraphDataVariantT;

template struct Node_t<Kmer<64>::Type >; 
template struct Edge_t<Node_t<Kmer<64>::Type > >; 
template class GraphUnitigsTemplate<64>; 
template class UnitigsConstructionAlgorithm<64>; 


template class Simplifications <GraphUnitigsTemplate<64>, NodeGU, EdgeGU >; 




/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
