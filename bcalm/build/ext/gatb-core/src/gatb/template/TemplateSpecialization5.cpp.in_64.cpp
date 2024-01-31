
// since we didn't define the functions in a .h file, that trick removes linker errors,
// see http://www.parashift.com/c++-faq-lite/separate-template-class-defn-from-decl.html

#include <gatb/debruijn/impl/BranchingAlgorithm.cpp>
#include <gatb/debruijn/impl/IterativeExtensions.cpp>
#include <gatb/debruijn/impl/Graph.cpp>
#include <gatb/debruijn/impl/GraphUnitigs.hpp>

/********************************************************************************/
namespace gatb { namespace core { namespace debruijn { namespace impl  {
/********************************************************************************/

template class BranchingAlgorithm   <64, Node, Edge, Graph>;

typedef boost::variant<GraphData<64> > GraphDataVariantT; // same as GraphDataVariantFast in Graph.hpp

template class BranchingAlgorithm <64,Node_t<Kmer<64>::Type>,Edge_t<Node_t<Kmer<64>::Type > >, GraphTemplate<Node_t<Kmer<64>::Type>,Edge_t<Node_t<Kmer<64>::Type > >,GraphDataVariantT>>;
//template class BranchingAlgorithm <64, NodeGU, EdgeGU, GraphUnitigsTemplate<64>>; // same as below

template class IterativeExtensions  <64, Node, Edge, Graph>;

template class IterativeExtensions  <64, Node_t<Kmer<64>::Type>,Edge_t<Node_t<Kmer<64>::Type > >, GraphTemplate<Node_t<Kmer<64>::Type>,Edge_t<Node_t<Kmer<64>::Type > >,GraphDataVariantT>>;
//template class IterativeExtensions  <64, NodeGU, EdgeGU, GraphUnitigsTemplate<64>>; // IterativeExtensinos isn't ready for GraphU, because it asks for Node::Value then queries node.kmer

/********************************************************************************/
} } } } /* end of namespaces. */
/********************************************************************************/
