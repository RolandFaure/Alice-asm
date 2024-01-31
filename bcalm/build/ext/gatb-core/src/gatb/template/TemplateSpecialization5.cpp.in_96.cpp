
// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

#include <gatb/debruijn/impl/BranchingAlgorithm.cpp>
#include <gatb/debruijn/impl/IterativeExtensions.cpp>
#include <gatb/debruijn/impl/Graph.cpp>
#include <gatb/debruijn/impl/GraphUnitigs.hpp>

/********************************************************************************/
namespace gatb { namespace core { namespace debruijn { namespace impl  {
/********************************************************************************/

template class BranchingAlgorithm   <96, Node, Edge, Graph>;

typedef boost::variant<GraphData<96> > GraphDataVariantT; // same as GraphDataVariantFast in Graph.hpp

template class BranchingAlgorithm <96,Node_t<Kmer<96>::Type>,Edge_t<Node_t<Kmer<96>::Type > >, GraphTemplate<Node_t<Kmer<96>::Type>,Edge_t<Node_t<Kmer<96>::Type > >,GraphDataVariantT>>;
//template class BranchingAlgorithm <96, NodeGU, EdgeGU, GraphUnitigsTemplate<96>>; // same as below

template class IterativeExtensions  <96, Node, Edge, Graph>;

template class IterativeExtensions  <96, Node_t<Kmer<96>::Type>,Edge_t<Node_t<Kmer<96>::Type > >, GraphTemplate<Node_t<Kmer<96>::Type>,Edge_t<Node_t<Kmer<96>::Type > >,GraphDataVariantT>>;
//template class IterativeExtensions  <96, NodeGU, EdgeGU, GraphUnitigsTemplate<96>>; // IterativeExtensinos isn't ready for GraphU, because it asks for Node::Value then queries node.kmer

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
