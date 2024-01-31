#include <gatb/debruijn/impl/Simplifications.cpp>
#include <gatb/debruijn/impl/UnitigsConstructionAlgorithm.cpp>
#include <gatb/debruijn/impl/GraphUnitigs.cpp>

using namespace gatb::core::kmer;
using namespace gatb::core::kmer::impl;

/********************************************************************************/
namespace gatb { namespace core { namespace debruijn { namespace impl  {
/********************************************************************************/

typedef boost::variant<GraphData<128> > GraphDataVariantT;

template struct Node_t<Kmer<128>::Type >; 
template struct Edge_t<Node_t<Kmer<128>::Type > >; 
template class GraphUnitigsTemplate<128>; 
template class UnitigsConstructionAlgorithm<128>; 


template class Simplifications <GraphUnitigsTemplate<128>, NodeGU, EdgeGU >; 




/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
