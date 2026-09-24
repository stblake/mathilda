/* test_graph_algos.c - the "algos" stream of the graph subsystem: flows and
 * cuts, matchings and covers, independent sets and cliques, Hamiltonian
 * cycles/paths, isomorphism and planarity (src/graph/galg_*.c).
 *
 * Expected values were cross-checked against Mathematica 15 (identical output
 * wherever the answer is unique; where Mathematica's choice among equally
 * optimal answers differs, the test pins Mathilda's deterministic choice and a
 * separate assertion checks validity). The randomized differential test
 * against Mathematica lives in benchmarks/95-graph-algorithms/diff_mathematica.py.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "print.h"
#include "graph.h"
#include "graph_algos.h"
#include "test_utils.h"
#include <stdlib.h>

#define G4 "Graph[{1<->2,2<->3,3<->4,4<->1,1<->3}]"
#define GW "Graph[{1,2,3,4},{1<->2,2<->3,3<->4,4<->1,1<->3},EdgeWeight->{5,5,10,10,5}]"
#define HD "Graph[{1->2,1->3,2->3,2->4,3->4}]"

static void test_max_flow(void) {
    assert_eval_eq("FindMaximumFlow[" G4 ",1,3]", "3", 0);
    assert_eval_eq("FindMaximumFlow[" G4 ",1,3,\"FlowMatrix\"]",
                   "{{0, 1, 1, 1}, {0, 0, 1, 0}, {0, 0, 0, 0}, {0, 0, 1, 0}}", 0);
    assert_eval_eq("FindMaximumFlow[" G4 ",1,3,\"EdgeList\"]",
                   "{1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 3, 4 <-> 3}", 0);
    assert_eval_eq("FindMaximumFlow[" G4 ",1,3,\"FlowValue\"]", "3", 0);
    /* EdgeCapacity option; EdgeWeight is ignored by FindMaximumFlow */
    assert_eval_eq("FindMaximumFlow[" HD ",1,4,EdgeCapacity->{3,2,1,2,3}]", "5", 0);
    assert_eval_eq("FindMaximumFlow[Graph[{1,2,3,4},{1->2,1->3,2->3,2->4,3->4},"
                   "EdgeWeight->{3,2,1,2,3}],1,4]", "2", 0);
    assert_eval_eq("FindMaximumFlow[" HD ",1,4,\"FlowMatrix\",EdgeCapacity->{3,2,1,2,3}]",
                   "{{0, 3, 2, 0}, {0, 0, 1, 2}, {0, 0, 0, 3}, {0, 0, 0, 0}}", 0);
    assert_eval_eq("FindMaximumFlow[" HD ",1,4,\"EdgeList\",EdgeCapacity->{3,2,1,2,3}]",
                   "{1 -> 2, 1 -> 3, 2 -> 3, 2 -> 4, 3 -> 4}", 0);
    /* numeric kinds */
    assert_eval_eq("FindMaximumFlow[Graph[{1->2,2->3}],1,3,EdgeCapacity->{1.5,2}]", "1.5", 0);
    assert_eval_eq("FindMaximumFlow[Graph[{1->2,2->3}],1,3,EdgeCapacity->{3/2,2}]", "1.5", 0);
    assert_eval_eq("FindMaximumFlow[Graph[{1->2,2->3}],1,3,EdgeCapacity->{Infinity,2}]", "2", 0);
    assert_eval_eq("FindMaximumFlow[Graph[{1->2,2->3}],1,3,EdgeCapacity->{Infinity,Infinity}]",
                   "Infinity", 0);
    assert_eval_eq("Head[FindMaximumFlow[Graph[{1->2,2->3}],1,3,EdgeCapacity->{x,2}]]",
                   "FindMaximumFlow", 0);
    /* terminals */
    assert_eval_eq("FindMaximumFlow[" G4 ",1,1]", "0", 0);
    assert_eval_eq("FindMaximumFlow[" G4 ",{1,2},{3,4}]", "3", 0);
    assert_eval_eq("FindMaximumFlow[Graph[{1->2,2->3,3->4,1->3}],{1,2},4]", "1", 0);
    assert_eval_eq("FindMaximumFlow[Graph[{1<->2,2->3}],3,1]", "0", 0);
    assert_eval_eq("Head[FindMaximumFlow[" G4 ",1,5]]", "FindMaximumFlow", 0);
    /* vertex capacities */
    assert_eval_eq("FindMaximumFlow[Graph[{1->2,2->3,3->4,1->3}],1,4,"
                   "VertexCapacity->{Infinity,1,1,Infinity}]", "1", 0);
    assert_eval_eq("FindMaximumFlow[Graph[{1->2,2->3,3->4,1->3}],1,4,"
                   "VertexCapacity->{1,1,1,1}]", "1", 0);
    /* packed flow matrix */
    assert_eval_eq("NDArrayQ[FindMaximumFlow[CycleGraph[4],1,3,\"FlowMatrix\"]]", "True", 0);
}

static void test_cuts(void) {
    assert_eval_eq("FindMinimumCut[" G4 "]", "{2, {{2}, {1, 3, 4}}}", 0);
    assert_eval_eq("FindMinimumCut[" GW "]", "{10, {{2}, {1, 3, 4}}}", 0);
    assert_eval_eq("FindMinimumCut[Graph[{1,2,3},{1<->2}]]", "{0, {{3}, {1, 2}}}", 0);
    assert_eval_eq("FindMinimumCut[Graph[{1->2,2->3,3->1,1->3}]]", "{1, {{1, 3}, {2}}}", 0);
    assert_eval_eq("FindMinimumCut[Graph[{1,2,3},{1->2,2->3,3->1},EdgeWeight->{2,3,4}]]",
                   "{2, {{1}, {2, 3}}}", 0);
    assert_eval_eq("Head[FindMinimumCut[Graph[{1},{}]]]", "FindMinimumCut", 0);
    assert_eval_eq("FindEdgeCut[" GW "]", "{1 <-> 2, 2 <-> 3}", 0);
    assert_eval_eq("FindEdgeCut[" G4 ",1,3]", "{1 <-> 2, 4 <-> 1, 1 <-> 3}", 0);
    assert_eval_eq("FindEdgeCut[" GW ",1,3]", "{1 <-> 2, 4 <-> 1, 1 <-> 3}", 0);
    assert_eval_eq("FindEdgeCut[Graph[{1->2,2->3,3->4}],1,4]", "{1 -> 2}", 0);
    assert_eval_eq("FindEdgeCut[Graph[{1->2,2->3,3->4}],4,1]", "{}", 0);
    assert_eval_eq("FindEdgeCut[Graph[{1->2,2->3,3->4}]]", "{}", 0);
    assert_eval_eq("FindEdgeCut[Graph[{1},{}]]", "{}", 0);
    assert_eval_eq("EdgeConnectivity[" G4 "]", "2", 0);
    assert_eval_eq("EdgeConnectivity[" G4 ",1,3]", "3", 0);
    assert_eval_eq("EdgeConnectivity[" GW "]", "10", 0);
    assert_eval_eq("EdgeConnectivity[" GW ",1,3]", "20", 0);
    assert_eval_eq("EdgeConnectivity[Graph[{1,2,3,4},{1<->2,2<->3,3<->4,4<->1,1<->3},"
                   "EdgeWeight->{5,5,10,10,5/2}]]", "10.0", 0);
    assert_eval_eq("EdgeConnectivity[" HD "]", "0", 0);
    assert_eval_eq("EdgeConnectivity[Graph[{1->2,2->3,3->1,1->3}]]", "1", 0);
    assert_eval_eq("EdgeConnectivity[Graph[{1<->2,2->3,3<->1}]]", "1", 0);
    assert_eval_eq("Head[EdgeConnectivity[Graph[{1},{}]]]", "EdgeConnectivity", 0);
    /* vertex cuts: underlying undirected graph, Mathematica's conventions */
    assert_eval_eq("FindVertexCut[" G4 "]", "{1, 3}", 0);
    assert_eval_eq("FindVertexCut[" G4 ",2,4]", "{1, 3}", 0);
    assert_eval_eq("FindVertexCut[" G4 ",1,3]", "{}", 0);
    assert_eval_eq("FindVertexCut[CompleteGraph[4]]", "{1, 2, 3}", 0);
    assert_eval_eq("FindVertexCut[CycleGraph[3]]", "{1, 2}", 0);
    assert_eval_eq("FindVertexCut[Graph[{1,2,3},{1<->2}]]", "{}", 0);
    assert_eval_eq("FindVertexCut[PathGraph[{1,2,3,4}],1,4]", "{3}", 0);
    assert_eval_eq("FindVertexCut[Graph[{1->2,2->3,3->4,4->1}]]", "{2, 4}", 0);
    assert_eval_eq("FindVertexCut[Graph[{1->2,2->3,3->1,1->4,4->2,3->4}]]", "{}", 0);
    assert_eval_eq("FindVertexCut[Graph[{1->2,2->1}]]", "{1}", 0);
}

static void test_matching_covers(void) {
    assert_eval_eq("FindIndependentEdgeSet[" G4 "]", "{1 <-> 2, 3 <-> 4}", 0);
    assert_eval_eq("FindIndependentEdgeSet[Graph[{1->2,2->3,3->4}]]", "{1 -> 2, 3 -> 4}", 0);
    assert_eval_eq("FindIndependentEdgeSet[Graph[{1,2,3},{}]]", "{}", 0);
    assert_eval_eq("Length[FindIndependentEdgeSet[CycleGraph[5]]]", "2", 0);
    /* Petersen graph has a perfect matching (non-bipartite: blossoms) */
    assert_eval_eq("Length[FindIndependentEdgeSet[Graph[Range[10],{1<->2,2<->3,3<->4,4<->5,"
                   "5<->1,1<->6,2<->7,3<->8,4<->9,5<->10,6<->8,8<->10,10<->7,7<->9,9<->6}]]]", "5", 0);
    assert_eval_eq("FindEdgeCover[" G4 "]", "{1 <-> 2, 3 <-> 4}", 0);
    assert_eval_eq("FindEdgeCover[Graph[{1,2,3},{1<->2}]]", "{}", 0);
    assert_eval_eq("Length[FindEdgeCover[CycleGraph[5]]]", "3", 0);
    assert_eval_eq("FindEdgeCover[Graph[{1->2,2->3,3->4}]]", "{1 -> 2, 3 -> 4}", 0);
    assert_eval_eq("EdgeCoverQ[CycleGraph[7], FindEdgeCover[CycleGraph[7]]]", "True", 0);
    assert_eval_eq("FindVertexCover[" G4 "]", "{1, 3}", 0);
    assert_eval_eq("FindVertexCover[CycleGraph[5]]", "{1, 2, 4}", 0);
    assert_eval_eq("FindVertexCover[Graph[{1->2,2->3,3->4}]]", "{2, 4}", 0);
    assert_eval_eq("FindIndependentVertexSet[" G4 "]", "{{2, 4}}", 0);
    assert_eval_eq("FindIndependentVertexSet[StarGraph[5]]", "{{2, 3, 4, 5}}", 0);
    assert_eval_eq("FindIndependentVertexSet[CycleGraph[5],2]", "{{1, 3}}", 0);
    assert_eval_eq("FindIndependentVertexSet[CycleGraph[5],{2}]", "{{1, 3}}", 0);
    assert_eval_eq("FindIndependentVertexSet[CycleGraph[5],Infinity,All]",
                   "{{3, 5}, {2, 5}, {2, 4}, {1, 4}, {1, 3}}", 0);
    assert_eval_eq("FindIndependentVertexSet[CycleGraph[5],{2},3]",
                   "{{2, 4}, {1, 4}, {1, 3}}", 0);
    assert_eval_eq("FindIndependentVertexSet[StarGraph[5],Infinity,All]",
                   "{{2, 3, 4, 5}, {1}}", 0);
    assert_eval_eq("FindIndependentVertexSet[Graph[{1->2,2->3,3->4}]]", "{{1, 3}}", 0);
    /* exactness on a folding-heavy sparse graph: alpha(C_101) = 50 */
    assert_eval_eq("Length[First[FindIndependentVertexSet[CycleGraph[101]]]]", "50", 0);
    assert_eval_eq("Length[FindVertexCover[CycleGraph[101]]]", "51", 0);
    /* Petersen: alpha = 4 */
    assert_eval_eq("Length[First[FindIndependentVertexSet[Graph[Range[10],{1<->2,2<->3,3<->4,"
                   "4<->5,5<->1,1<->6,2<->7,3<->8,4<->9,5<->10,6<->8,8<->10,10<->7,7<->9,9<->6}]]]]",
                   "4", 0);
}

static void test_set_predicates(void) {
    assert_eval_eq("{IndependentVertexSetQ[" G4 ",{2,4}], IndependentVertexSetQ[" G4 ",{1,4}], "
                   "IndependentVertexSetQ[" G4 ",{}], IndependentVertexSetQ[" G4 ",{5}], "
                   "IndependentVertexSetQ[" G4 ",{2,2}]}",
                   "{True, False, True, False, True}", 0);
    assert_eval_eq("{VertexCoverQ[" G4 ",{1,3}], VertexCoverQ[" G4 ",{1}], VertexCoverQ[" G4 ",{1,3,7}]}",
                   "{True, False, False}", 0);
    assert_eval_eq("{IndependentEdgeSetQ[" G4 ",{1<->2,3<->4}], IndependentEdgeSetQ[" G4 ",{2<->1,3<->4}], "
                   "IndependentEdgeSetQ[" G4 ",{1<->2,1<->3}], IndependentEdgeSetQ[" G4 ",{1<->5}], "
                   "IndependentEdgeSetQ[" G4 ",{}]}",
                   "{True, True, False, False, True}", 0);
    assert_eval_eq("{EdgeCoverQ[" G4 ",{1<->2,3<->4}], EdgeCoverQ[" G4 ",{1<->2}], "
                   "EdgeCoverQ[Graph[{1->2,3->4}],{1->2,3->4}], EdgeCoverQ[Graph[{1->2,3->4}],{2->1,3->4}], "
                   "IndependentEdgeSetQ[Graph[{1->2,3->4}],{2->1}]}",
                   "{True, False, True, False, False}", 0);
    assert_eval_eq("IndependentVertexSetQ[5,{1}]", "False", 0);
}

static void test_cliques(void) {
    assert_eval_eq("FindClique[" G4 "]", "{{1, 2, 3}}", 0);
    assert_eval_eq("FindClique[" G4 ",2]", "{}", 0);
    assert_eval_eq("FindClique[" G4 ",{2}]", "{}", 0);
    assert_eval_eq("FindClique[" G4 ",{5}]", "{}", 0);
    assert_eval_eq("FindClique[" G4 ",5]", "{{1, 2, 3}}", 0);
    assert_eval_eq("FindClique[" G4 ",Infinity,All]", "{{1, 3, 4}, {1, 2, 3}}", 0);
    assert_eval_eq("FindClique[" G4 ",3,All]", "{{1, 3, 4}, {1, 2, 3}}", 0);
    assert_eval_eq("FindClique[StarGraph[5],Infinity,All]", "{{1, 5}, {1, 4}, {1, 3}, {1, 2}}", 0);
    assert_eval_eq("FindClique[CycleGraph[5],Infinity,All]",
                   "{{4, 5}, {3, 4}, {2, 3}, {1, 5}, {1, 2}}", 0);
    assert_eval_eq("FindClique[Graph[{1->2,2->1,2->3}]]", "{{1, 2}}", 0);
    assert_eval_eq("FindClique[Graph[{1->2,2->3,3->1}]]", "{{1}}", 0);
    assert_eval_eq("FindClique[Graph[{1,2,3},{}]]", "{{1}}", 0);
    assert_eval_eq("FindClique[Graph[{},{}]]", "{}", 0);
    assert_eval_eq("FindClique[CompleteGraph[40]] === {Range[40]}", "True", 0);
    assert_eval_eq("FindKClique[" G4 ",2]", "{{1, 2, 3, 4}}", 0);
    assert_eval_eq("FindKClique[CycleGraph[8],2]", "{{1, 2, 3}}", 0);
}

static void test_hamiltonian(void) {
    assert_eval_eq("FindHamiltonianCycle[" G4 "]", "{{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}}", 0);
    assert_eval_eq("FindHamiltonianCycle[StarGraph[5]]", "{}", 0);
    assert_eval_eq("FindHamiltonianCycle[Graph[{1->2,2->3,3->1}]]", "{{1 -> 2, 2 -> 3, 3 -> 1}}", 0);
    assert_eval_eq("FindHamiltonianCycle[Graph[{1->2,2->1}]]", "{{1 -> 2, 2 -> 1}}", 0);
    assert_eval_eq("FindHamiltonianCycle[Graph[{1,2},{1<->2}]]", "{}", 0);
    assert_eval_eq("FindHamiltonianCycle[Graph[{1},{}]]", "{{}}", 0);
    assert_eval_eq("Length[FindHamiltonianCycle[CompleteGraph[4],All]]", "3", 0);
    assert_eval_eq("Length[FindHamiltonianCycle[CompleteGraph[6],All]]", "60", 0);
    assert_eval_eq("Length[FindHamiltonianCycle[CompleteGraph[6],7]]", "7", 0);
    assert_eval_eq("FindHamiltonianPath[StarGraph[5]]", "{}", 0);
    assert_eval_eq("FindHamiltonianPath[Graph[{1->2,2->3,3->1}]]", "{1, 2, 3}", 0);
    assert_eval_eq("FindHamiltonianPath[Graph[{1,2,3},{1<->2,2<->3}],3,1]", "{3, 2, 1}", 0);
    assert_eval_eq("FindHamiltonianPath[Graph[{1,2,3},{1<->2,2<->3}],1,2]", "{}", 0);
    assert_eval_eq("FindHamiltonianPath[Graph[{1},{}]]", "{}", 0);
    assert_eval_eq("{HamiltonianGraphQ[" G4 "], HamiltonianGraphQ[StarGraph[5]], "
                   "HamiltonianGraphQ[Graph[{1},{}]], HamiltonianGraphQ[Graph[{1,2},{1<->2}]], "
                   "HamiltonianGraphQ[Graph[{},{}]], HamiltonianGraphQ[7]}",
                   "{True, False, True, False, False, False}", 0);
    /* Petersen: no Hamiltonian cycle, but a Hamiltonian path */
    assert_eval_eq("HamiltonianGraphQ[Graph[Range[10],{1<->2,2<->3,3<->4,4<->5,5<->1,1<->6,2<->7,"
                   "3<->8,4<->9,5<->10,6<->8,8<->10,10<->7,7<->9,9<->6}]]", "False", 0);
    assert_eval_eq("Length[FindHamiltonianPath[Graph[Range[10],{1<->2,2<->3,3<->4,4<->5,5<->1,"
                   "1<->6,2<->7,3<->8,4<->9,5<->10,6<->8,8<->10,10<->7,7<->9,9<->6}]]]", "10", 0);
}

static void test_planarity(void) {
    assert_eval_eq("PlanarGraphQ[CompleteGraph[4]]", "True", 0);
    assert_eval_eq("PlanarGraphQ[CompleteGraph[5]]", "False", 0);
    assert_eval_eq("PlanarGraphQ[Graph[{1->2,2->1}]]", "True", 0);
    assert_eval_eq("PlanarGraphQ[Graph[{},{}]]", "True", 0);
    assert_eval_eq("PlanarGraphQ[5]", "False", 0);
    /* K3,3 */
    assert_eval_eq("PlanarGraphQ[Graph[{1<->4,1<->5,1<->6,2<->4,2<->5,2<->6,3<->4,3<->5,3<->6}]]",
                   "False", 0);
    assert_eval_eq("PlanarGraphQ[Graph[Range[10],{1<->2,2<->3,3<->4,4<->5,5<->1,1<->6,2<->7,"
                   "3<->8,4<->9,5<->10,6<->8,8<->10,10<->7,7<->9,9<->6}]]", "False", 0);
}

static void test_isomorphism(void) {
    assert_eval_eq("IsomorphicGraphQ[CycleGraph[5], Graph[{1<->3,3<->5,5<->2,2<->4,4<->1}]]", "True", 0);
    assert_eval_eq("IsomorphicGraphQ[CycleGraph[6], Graph[{1<->2,2<->3,3<->1,4<->5,5<->6,6<->4}]]", "False", 0);
    assert_eval_eq("IsomorphicGraphQ[PathGraph[Range[4]], PathGraph[{4,2,3,1}], PathGraph[{a,b,c,d}]]", "True", 0);
    assert_eval_eq("IsomorphicGraphQ[CycleGraph[3], 5]", "False", 0);
    assert_eval_eq("IsomorphicGraphQ[Graph[{1->2}], Graph[{1<->2}]]", "False", 0);
    assert_eval_eq("IsomorphicGraphQ[Graph[{1->2,2->3}], Graph[{3->1,2->3}]]", "True", 0);
    assert_eval_eq("IsomorphicGraphQ[Graph[{1->2,2->3}], Graph[{1->2,3->2}]]", "False", 0);
    assert_eval_eq("IsomorphicGraphQ[Graph[{1->2,2<->3}], Graph[{3<->1,2->3}]]", "True", 0);
    assert_eval_eq("IsomorphicGraphQ[Graph[{1->2,2<->3}], Graph[{3<->1,3->2}]]", "False", 0);
    assert_eval_eq("IsomorphicGraphQ[Graph[{},{}], Graph[{},{}]]", "True", 0);
    /* Petersen graph vs a relabeled copy */
    assert_eval_eq("IsomorphicGraphQ[Graph[Range[10],{1<->2,2<->3,3<->4,4<->5,5<->1,1<->6,2<->7,"
                   "3<->8,4<->9,5<->10,6<->8,8<->10,10<->7,7<->9,9<->6}], Graph[Range[10],"
                   "{10<->9,9<->8,8<->7,7<->6,6<->10,10<->5,9<->4,8<->3,7<->2,6<->1,5<->3,3<->1,1<->4,4<->2,2<->5}]]",
                   "True", 0);

    assert_eval_eq("FindGraphIsomorphism[Graph[{3,1,2},{1<->2,2<->3}], Graph[{a,b,c},{a<->b,a<->c}], All]",
                   "{<|3 -> b, 1 -> c, 2 -> a|>, <|3 -> c, 1 -> b, 2 -> a|>}", 0);
    assert_eval_eq("FindGraphIsomorphism[CycleGraph[3], PathGraph[{1,2,3}]]", "{}", 0);
    assert_eval_eq("Length[FindGraphIsomorphism[CycleGraph[6], CycleGraph[6], All]]", "12", 0);
    assert_eval_eq("Length[FindGraphIsomorphism[CompleteGraph[5], CompleteGraph[5], All]]", "120", 0);
    assert_eval_eq("Length[FindGraphIsomorphism[CycleGraph[6], CycleGraph[6], 5]]", "5", 0);
    assert_eval_eq("With[{h = Graph[{1<->3,3<->5,5<->7,7<->2,2<->4,4<->6,6<->1}]}, With[{a = First[FindGraphIsomorphism[CycleGraph[7], h]]},"
                   " Sort[Sort /@ (EdgeList[CycleGraph[7]] /. Normal[a])] === Sort[Sort /@ EdgeList[h]]]]",
                   "True", 0);

    assert_eval_eq("CanonicalGraph[Graph[{a,b,c},{a<->b,b<->c}]] === CanonicalGraph[Graph[{x,y,z},{z<->x,y<->z}]]", "True", 0);
    assert_eval_eq("CanonicalGraph[CycleGraph[6]] === CanonicalGraph[Graph[{1<->2,2<->3,3<->1,4<->5,5<->6,6<->4}]]", "False", 0);
    assert_eval_eq("VertexList[CanonicalGraph[Graph[{a,b,c},{a->b,c->b}]]]", "{1, 2, 3}", 0);
    assert_eval_eq("CanonicalGraph[Graph[{a,b,c},{a->b,c->b}]] === CanonicalGraph[Graph[{a,b,c},{b->a,c->a}]]", "True", 0);

    assert_eval_eq("GraphAutomorphismGroup[Graph[{c,a,b},{c<->a}]]", "PermutationGroup[{Cycles[{{1, 2}}]}]", 0);
    assert_eval_eq("GraphAutomorphismGroup[Graph[{1,2},{1->2}]]", "PermutationGroup[{}]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_max_flow);
    TEST(test_cuts);
    TEST(test_matching_covers);
    TEST(test_set_predicates);
    TEST(test_cliques);
    TEST(test_hamiltonian);
    TEST(test_planarity);
    TEST(test_isomorphism);

    printf("All graph algos tests passed!\n");
    return 0;
}
