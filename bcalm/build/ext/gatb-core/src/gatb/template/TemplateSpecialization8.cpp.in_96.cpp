#include <gatb/debruijn/impl/Simplifications.cpp>
#include <gatb/debruijn/impl/UnitigsConstructionAlgorithm.cpp>
#include <gatb/debruijn/impl/GraphUnitigs.cpp>

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

/********************************************************************************/
namespace gatb { namespace core { namespace debruijn { namespace impl  {
/********************************************************************************/

typedef boost::variant<GraphData<96> > GraphDataVariantT;

template struct Node_t<Kmer<96>::Type >; 
template struct Edge_t<Node_t<Kmer<96>::Type > >; 
template class GraphUnitigsTemplate<96>; 
template class UnitigsConstructionAlgorithm<96>; 


template class Simplifications <GraphUnitigsTemplate<96>, NodeGU, EdgeGU >; 




/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
