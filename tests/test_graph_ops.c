/* test_graph_ops.c - graph editing, transforms, set operations, predicates,
 * Eulerian cycles, FindCycle and FindPath (src/graph/gops_*.c).
 *
 * Expected values were checked against Mathematica 15 wherever Mathilda
 * follows it exactly (orders included); the documented deviations (multigraph
 * results left unevaluated, FindCycle's choice of cycle for length-bounded
 * searches, directed LineGraph numbering) are tested as Mathilda's own
 * contract.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "print.h"
#include "graph.h"
#include "test_utils.h"
#include <stdlib.h>

#define G4 "Graph[{1,2,3,4},{1<->2,2<->3,3<->4}]"
#define W4 "Graph[{1,2,3,4},{1<->2,2<->3,3<->4},EdgeWeight->{5,6,7}]"
#define D3 "Graph[{1,2,3},{1->2,2->3}]"

static void setup(void) {
    assert_eval_eq("gv[h_] := {VertexList[h], EdgeList[h]}", "Null", 0);
    assert_eval_eq("gw[h_] := {VertexList[h], EdgeList[h], EdgeWeight[h]}", "Null", 0);
}

static void test_vertex_add_delete(void) {
    assert_eval_eq("gv[VertexAdd[" G4 ", 5]]", "{{1, 2, 3, 4, 5}, {1 <-> 2, 2 <-> 3, 3 <-> 4}}", 0);
    assert_eval_eq("gv[VertexAdd[" G4 ", {5, 1, 6, 5}]]",
                   "{{1, 2, 3, 4, 5, 6}, {1 <-> 2, 2 <-> 3, 3 <-> 4}}", 0);
    assert_eval_eq("gw[VertexAdd[" W4 ", 7]]",
                   "{{1, 2, 3, 4, 7}, {1 <-> 2, 2 <-> 3, 3 <-> 4}, {5, 6, 7}}", 0);
    assert_eval_eq("gv[VertexAdd[Graph[{1,2},{1<->2}], {{3,4}}]]", "{{1, 2, {3, 4}}, {1 <-> 2}}", 0);
    assert_eval_eq("gv[VertexAdd[Graph[{1,2},{1<->2}], 1]]", "{{1, 2}, {1 <-> 2}}", 0);
    /* weights stay aligned with the surviving edges */
    assert_eval_eq("gw[VertexDelete[" W4 ", 2]]", "{{1, 3, 4}, {3 <-> 4}, {7}}", 0);
    assert_eval_eq("gv[VertexDelete[Graph[{5,4,3,2,1},{4<->3,1<->2,5<->4,2<->3}], 1]]",
                   "{{5, 4, 3, 2}, {4 <-> 3, 5 <-> 4, 2 <-> 3}}", 0);
    assert_eval_eq("gv[VertexDelete[Graph[{1,2,3},{1<->2,2<->3}], {2,2}]]", "{{1, 3}, {}}", 0);
    assert_eval_eq("gv[VertexDelete[" D3 ", _?EvenQ]]", "{{1, 3}, {}}", 0);
    /* a listed non-vertex leaves the call unevaluated */
    assert_eval_eq("Head[VertexDelete[" G4 ", {2, 9}]]", "VertexDelete", 0);
    assert_eval_eq("Head[VertexDelete[" G4 ", 9]]", "VertexDelete", 0);
}

static void test_edge_add_delete(void) {
    assert_eval_eq("gv[EdgeAdd[" G4 ", 1<->4]]",
                   "{{1, 2, 3, 4}, {1 <-> 2, 2 <-> 3, 3 <-> 4, 1 <-> 4}}", 0);
    /* u->v sugar takes the graph's kind; new endpoints become vertices */
    assert_eval_eq("gv[EdgeAdd[" G4 ", {1->5, 6<->7}]]",
                   "{{1, 2, 3, 4, 5, 6, 7}, {1 <-> 2, 2 <-> 3, 3 <-> 4, 1 <-> 5, 6 <-> 7}}", 0);
    assert_eval_eq("gv[EdgeAdd[Graph[{1,2},{}], 1->2]]", "{{1, 2}, {1 <-> 2}}", 0);
    assert_eval_eq("gv[EdgeAdd[Graph[{1,2},{}], DirectedEdge[1,2]]]", "{{1, 2}, {1 -> 2}}", 0);
    assert_eval_eq("gv[EdgeAdd[Graph[{1,2,3},{1<->2, 2->3}], 3->1]]",
                   "{{1, 2, 3}, {1 <-> 2, 2 -> 3, 3 -> 1}}", 0);
    assert_eval_eq("gv[EdgeAdd[" D3 ", 2->1]]", "{{1, 2, 3}, {1 -> 2, 2 -> 3, 2 -> 1}}", 0);
    assert_eval_eq("gw[EdgeAdd[" W4 ", 1<->4]]",
                   "{{1, 2, 3, 4}, {1 <-> 2, 2 <-> 3, 3 <-> 4, 1 <-> 4}, {5, 6, 7, 1}}", 0);
    /* multigraph / self-loop results: unevaluated */
    assert_eval_eq("Head[EdgeAdd[" G4 ", 2<->1]]", "EdgeAdd", 0);
    assert_eval_eq("Head[EdgeAdd[" G4 ", 1<->1]]", "EdgeAdd", 0);
    assert_eval_eq("Head[EdgeAdd[Graph[{1,2},{1<->2}], {2<->3, 3<->2}]]", "EdgeAdd", 0);

    assert_eval_eq("gw[EdgeDelete[" W4 ", 2<->3]]", "{{1, 2, 3, 4}, {1 <-> 2, 3 <-> 4}, {5, 7}}", 0);
    assert_eval_eq("gv[EdgeDelete[" G4 ", 3<->2]]", "{{1, 2, 3, 4}, {1 <-> 2, 3 <-> 4}}", 0);
    assert_eval_eq("gv[EdgeDelete[Graph[{1,2,3},{1<->2,2<->3}], {2<->1,1<->2}]]",
                   "{{1, 2, 3}, {2 <-> 3}}", 0);
    assert_eval_eq("gv[EdgeDelete[" D3 ", _DirectedEdge]]", "{{1, 2, 3}, {}}", 0);
    assert_eval_eq("Head[EdgeDelete[" G4 ", {1<->2, 9<->10}]]", "EdgeDelete", 0);
    assert_eval_eq("Head[EdgeDelete[" D3 ", 2->1]]", "EdgeDelete", 0);
    assert_eval_eq("Head[EdgeDelete[" G4 ", 1->2]]", "EdgeDelete", 0);
}

static void test_subgraph_neighborhood(void) {
    assert_eval_eq("gw[Subgraph[" W4 ", {3, 2, 4}]]", "{{3, 2, 4}, {2 <-> 3, 3 <-> 4}, {6, 7}}", 0);
    assert_eval_eq("gv[Subgraph[Graph[{1,2,3,4},{3<->4,1<->2,2<->3}], {3,2,2,4}]]",
                   "{{3, 2, 4}, {2 <-> 3, 3 <-> 4}}", 0);
    assert_eval_eq("gv[Subgraph[Graph[{1,2,3,4,5},{4->3,3->4,1<->2,5<->4,2<->3}], {5,4,3,2}]]",
                   "{{5, 4, 3, 2}, {5 <-> 4, 3 -> 4, 4 -> 3, 2 <-> 3}}", 0);
    assert_eval_eq("gv[Subgraph[Graph[{1,2,3,4,5},{4->3,1<->2,5<->4,2<->3}], {1,2,3,4,5}]]",
                   "{{1, 2, 3, 4, 5}, {1 <-> 2, 2 <-> 3, 4 -> 3, 5 <-> 4}}", 0);
    assert_eval_eq("gv[Subgraph[" G4 ", {1, 9}]]", "{{1}, {}}", 0);
    assert_eval_eq("gv[Subgraph[" G4 ", {}]]", "{{}, {}}", 0);

    assert_eval_eq("gv[NeighborhoodGraph[" G4 ", 2]]", "{{2, 1, 3}, {1 <-> 2, 2 <-> 3}}", 0);
    assert_eval_eq("gv[NeighborhoodGraph[" G4 ", {1,4}]]",
                   "{{1, 4, 2, 3}, {1 <-> 2, 2 <-> 3, 3 <-> 4}}", 0);
    assert_eval_eq("gv[NeighborhoodGraph[" D3 ", 2]]", "{{2, 1, 3}, {1 -> 2, 2 -> 3}}", 0);
    assert_eval_eq("gv[NeighborhoodGraph[Graph[{1,2,3,4,5,6,7},{1<->5,1<->3,4<->1,3<->6,5<->7,2<->6}], 1, 2]]",
                   "{{1, 3, 4, 5, 6, 7}, {1 <-> 3, 4 <-> 1, 1 <-> 5, 3 <-> 6, 5 <-> 7}}", 0);
    assert_eval_eq("gv[NeighborhoodGraph[Graph[{1,2,3,4,5,6,7},{1<->5,1<->3,4<->1,3<->6,5<->7,2<->6}], {6,5}]]",
                   "{{6, 5, 2, 3, 1, 7}, {2 <-> 6, 3 <-> 6, 1 <-> 5, 1 <-> 3, 5 <-> 7}}", 0);
    assert_eval_eq("gv[NeighborhoodGraph[Graph[{1,2,3,4,5},{1->2,3->1,4->3,2->5}], 1, 2]]",
                   "{{1, 2, 3, 4, 5}, {1 -> 2, 3 -> 1, 4 -> 3, 2 -> 5}}", 0);
    assert_eval_eq("gv[NeighborhoodGraph[" G4 ", 1, 0]]", "{{1}, {}}", 0);
    assert_eval_eq("gv[NeighborhoodGraph[" G4 ", 7]]", "{{}, {}}", 0);
    assert_eval_eq("VertexCount[NeighborhoodGraph[PathGraph[Range[50]], 1, Infinity]]", "50", 0);
}

static void test_replace_rules_index(void) {
    assert_eval_eq("gw[VertexReplace[" W4 ", {1 -> a, 2 -> b}]]",
                   "{{a, b, 3, 4}, {a <-> b, b <-> 3, 3 <-> 4}, {5, 6, 7}}", 0);
    assert_eval_eq("gv[VertexReplace[Graph[{1,2,3},{1<->2,2<->3}], {x_ :> x^2, 3 -> 5}]]",
                   "{{1, 4, 9}, {1 <-> 4, 4 <-> 9}}", 0);
    assert_eval_eq("gv[VertexReplace[Graph[{1,2,3},{1<->2,2<->3}], {3->1, 1->3}]]",
                   "{{3, 2, 1}, {3 <-> 2, 2 <-> 1}}", 0);
    assert_eval_eq("gv[VertexReplace[" D3 ", {1 -> 3}]]", "{{3, 2}, {3 -> 2, 2 -> 3}}", 0);
    assert_eval_eq("Head[VertexReplace[" G4 ", 1 -> 2]]", "VertexReplace", 0);  /* loop */

    assert_eval_eq("EdgeRules[" G4 "]", "{1 -> 2, 2 -> 3, 3 -> 4}", 0);
    assert_eval_eq("EdgeRules[Graph[{1,2},{1<->2},EdgeWeight->{3}]]", "{1 -> 2}", 0);
    assert_eval_eq("{VertexIndex[" G4 ", 3], VertexIndex[" G4 ", {1, 3}]}", "{3, {1, 3}}", 0);
    assert_eval_eq("{EdgeIndex[" G4 ", 3<->2], EdgeIndex[Graph[{1,2},{1<->2}], {1<->2}], "
                   "EdgeIndex[Graph[{1,2,3},{1<->2}], 1->2]}", "{2, {1}, 1}", 0);
    assert_eval_eq("Head[VertexIndex[" G4 ", 9]]", "VertexIndex", 0);
    assert_eval_eq("Head[EdgeIndex[" D3 ", 2->1]]", "EdgeIndex", 0);

    assert_eval_eq("gv[IndexGraph[Graph[{a,b,c},{a->b,b->c}]]]", "{{1, 2, 3}, {1 -> 2, 2 -> 3}}", 0);
    assert_eval_eq("gv[IndexGraph[Graph[{a,b,c},{a->b,b->c}], 10]]",
                   "{{10, 11, 12}, {10 -> 11, 11 -> 12}}", 0);
    assert_eval_eq("gw[IndexGraph[" W4 ", 0]]",
                   "{{0, 1, 2, 3}, {0 <-> 1, 1 <-> 2, 2 <-> 3}, {5, 6, 7}}", 0);
    assert_eval_eq("Head[IndexGraph[Graph[{a,b},{a<->b}], x]]", "IndexGraph", 0);
}

static void test_transforms(void) {
    assert_eval_eq("gv[GraphComplement[" G4 "]]", "{{1, 2, 3, 4}, {1 <-> 3, 1 <-> 4, 2 <-> 4}}", 0);
    assert_eval_eq("gv[GraphComplement[" D3 "]]",
                   "{{1, 2, 3}, {1 -> 3, 2 -> 1, 3 -> 1, 3 -> 2}}", 0);
    assert_eval_eq("gv[GraphComplement[Graph[{1,2,3},{1->2, 2<->3}]]]",
                   "{{1, 2, 3}, {1 -> 3, 2 -> 1, 3 -> 1}}", 0);
    assert_eval_eq("gv[GraphComplement[Graph[{3,1,2},{3<->1}]]]", "{{3, 1, 2}, {3 <-> 2, 1 <-> 2}}", 0);
    assert_eval_eq("WeightedGraphQ[GraphComplement[" W4 "]]", "False", 0);

    assert_eval_eq("gw[ReverseGraph[Graph[{1,2,3},{2<->1, 3->2}, EdgeWeight->{4,5}]]]",
                   "{{1, 2, 3}, {2 <-> 1, 2 -> 3}, {4, 5}}", 0);
    assert_eval_eq("gv[UndirectedGraph[Graph[{1,2,3},{1->2,2->1,2->3}]]]",
                   "{{1, 2, 3}, {1 <-> 2, 2 <-> 3}}", 0);
    assert_eval_eq("gw[UndirectedGraph[Graph[{1,2,3},{2->1,1<->2,2->3}, EdgeWeight->{1,10,100}]]]",
                   "{{1, 2, 3}, {1 <-> 2, 2 <-> 3}, {11, 100}}", 0);
    assert_eval_eq("gv[UndirectedGraph[Graph[{3,2,1},{1->2, 3->1, 2->1}]]]",
                   "{{3, 2, 1}, {3 <-> 1, 2 <-> 1}}", 0);
    assert_eval_eq("gv[UndirectedGraph[Graph[{3,2,1},{3<->2, 1<->2}]]]",
                   "{{3, 2, 1}, {3 <-> 2, 1 <-> 2}}", 0);

    assert_eval_eq("gw[DirectedGraph[" W4 "]]",
                   "{{1, 2, 3, 4}, {1 -> 2, 2 -> 1, 2 -> 3, 3 -> 2, 3 -> 4, 4 -> 3}, {5, 5, 6, 6, 7, 7}}", 0);
    assert_eval_eq("gv[DirectedGraph[Graph[{1,2,3},{2<->1, 3->2}]]]",
                   "{{1, 2, 3}, {2 -> 1, 1 -> 2, 3 -> 2}}", 0);
    assert_eval_eq("gw[DirectedGraph[Graph[{4,2,1,3},{1<->2,3<->4,2<->3,1<->4,1<->3}, "
                   "EdgeWeight->{1,2,3,4,5}], \"Acyclic\"]]",
                   "{{4, 2, 1, 3}, {4 -> 1, 4 -> 3, 2 -> 1, 2 -> 3, 1 -> 3}, {4, 2, 1, 3, 5}}", 0);
    assert_eval_eq("gv[DirectedGraph[Graph[{4,2,1,3},{1->2, 3<->4, 2<->3}], \"Acyclic\"]]",
                   "{{4, 2, 1, 3}, {1 -> 2, 4 -> 3, 2 -> 3}}", 0);
    assert_eval_eq("AcyclicGraphQ[DirectedGraph[CompleteGraph[6], \"Acyclic\"]]", "True", 0);

    assert_eval_eq("gv[LineGraph[Graph[{1,2,3,4},{1<->2,2<->3,3<->4,4<->1,1<->3}]]]",
                   "{{1, 2, 3, 4, 5}, {2 <-> 1, 3 <-> 2, 4 <-> 3, 4 <-> 1, 5 <-> 1, 5 <-> 4, 5 <-> 2, 5 <-> 3}}", 0);
    assert_eval_eq("gv[LineGraph[Graph[{1,2,3},{1->2,2->3,3->1,1->3}]]]",
                   "{{1, 2, 3, 4}, {1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4, 4 -> 3}}", 0);
    assert_eval_eq("Head[LineGraph[Graph[{1,2,3},{1->2,2<->3}]]]", "LineGraph", 0);
}

static void test_set_operations(void) {
    assert_eval_eq("gv[GraphUnion[" D3 ", " G4 "]]",
                   "{{1, 2, 3, 4}, {1 -> 2, 2 -> 3, 1 <-> 2, 2 <-> 3, 3 <-> 4}}", 0);
    assert_eval_eq("gv[GraphUnion[Graph[{b<->a}], Graph[{c<->a}]]]", "{{a, b, c}, {a <-> b, a <-> c}}", 0);
    assert_eval_eq("gv[GraphUnion[Graph[{2<->1, 3->4}], Graph[{5<->1}]]]",
                   "{{1, 2, 3, 4, 5}, {3 -> 4, 2 <-> 1, 5 <-> 1}}", 0);
    assert_eval_eq("gv[GraphUnion[Graph[{4<->3, 2<->1}], Graph[{0<->5}]]]",
                   "{{0, 1, 2, 3, 4, 5}, {3 <-> 4, 1 <-> 2, 0 <-> 5}}", 0);
    assert_eval_eq("gv[GraphUnion[Graph[{4->3, 1->2}], Graph[{0->5}]]]",
                   "{{0, 1, 2, 3, 4, 5}, {4 -> 3, 1 -> 2, 0 -> 5}}", 0);
    assert_eval_eq("gv[GraphUnion[Graph[{1,2},{1<->2}], Graph[{2,1},{1->2}], Graph[{3,1},{3->1}]]]",
                   "{{1, 2, 3}, {1 -> 2, 3 -> 1, 1 <-> 2}}", 0);
    assert_eval_eq("gv[GraphUnion[Graph[{4,3,2,1},{4<->3, 2<->1}]]]",
                   "{{4, 3, 2, 1}, {4 <-> 3, 2 <-> 1}}", 0);
    assert_eval_eq("VertexList[GraphUnion[Graph[{b,1,\"s\",2.5,x[1]},{}], Graph[{a, 3/2}, {}]]]",
                   "{1, 3/2, 2.5, \"s\", a, b, x[1]}", 0);
    assert_eval_eq("WeightedGraphQ[GraphUnion[" W4 ", " W4 "]]", "False", 0);

    assert_eval_eq("gv[GraphIntersection[Graph[{4<->3, 1<->2}], Graph[{2<->1, 3<->4}]]]",
                   "{{1, 2, 3, 4}, {1 <-> 2, 4 <-> 3}}", 0);
    assert_eval_eq("gv[GraphIntersection[Graph[{4<->3, 1<->2, 5->6}], Graph[{2<->1, 3<->4, 5->6}]]]",
                   "{{1, 2, 3, 4, 5, 6}, {5 -> 6, 1 <-> 2, 4 <-> 3}}", 0);
    assert_eval_eq("gv[GraphIntersection[" G4 ", Graph[{5,3,2,1},{3<->2, 1<->5}]]]",
                   "{{1, 2, 3, 4, 5}, {2 <-> 3}}", 0);
    assert_eval_eq("gv[GraphDifference[Graph[{4<->3, 2<->1}], Graph[{5<->6}]]]",
                   "{{1, 2, 3, 4, 5, 6}, {2 <-> 1, 4 <-> 3}}", 0);
    assert_eval_eq("gv[GraphDifference[" G4 ", Graph[{3,2,9},{3<->2}]]]",
                   "{{1, 2, 3, 4, 9}, {1 <-> 2, 3 <-> 4}}", 0);
    assert_eval_eq("gv[GraphDisjointUnion[Graph[{a,b},{a->b}], Graph[{a,c},{a->c}]]]",
                   "{{1, 2, 3, 4}, {1 -> 2, 3 -> 4}}", 0);
    assert_eval_eq("gv[GraphDisjointUnion[Graph[{a->b}], Graph[{c<->d}], Graph[{x},{}]]]",
                   "{{1, 2, 3, 4, 5}, {1 -> 2, 3 <-> 4}}", 0);
    assert_eval_eq("gv[GraphDisjointUnion[Graph[{a,b},{a->b}]]]", "{{a, b}, {a -> b}}", 0);
    assert_eval_eq("Head[GraphDifference[" G4 "]]", "GraphDifference", 0);
}

static void test_predicates(void) {
    assert_eval_eq("{SimpleGraphQ[" G4 "], LoopFreeGraphQ[" G4 "], "
                   "MixedGraphQ[Graph[{1,2,3},{1->2, 2<->3}]], MixedGraphQ[" G4 "], "
                   "WeightedGraphQ[" W4 "], WeightedGraphQ[" G4 "], EdgeWeightedGraphQ[" W4 "], "
                   "SimpleGraphQ[5], MixedGraphQ[Graph[{1},{}]]}",
                   "{True, True, True, False, True, False, True, False, False}", 0);
    assert_eval_eq("{PathGraphQ[" G4 "], PathGraphQ[" D3 "], PathGraphQ[Graph[{1,2,3},{1->2,3->2}]], "
                   "PathGraphQ[Graph[{1},{}]], PathGraphQ[Graph[{},{}]], PathGraphQ[CycleGraph[3]], "
                   "PathGraphQ[Graph[{1,2,3,4},{1<->2,3<->4}]], PathGraphQ[Graph[{1,2,3},{1->2,2<->3}]]}",
                   "{True, True, False, True, False, True, False, False}", 0);
    assert_eval_eq("{PathGraphQ[Graph[{1->2,2->3,3->1}]], PathGraphQ[Graph[{1->2,2->1,2->3}]], "
                   "PathGraphQ[Graph[{2,1,3},{1<->2,2<->3}]], PathGraphQ[Graph[{1,2},{1->2,2->1}]]}",
                   "{True, False, True, True}", 0);
    assert_eval_eq("{EulerianGraphQ[CycleGraph[4]], EulerianGraphQ[" G4 "], "
                   "EulerianGraphQ[Graph[{1,2,3},{1->2,2->3,3->1}]], "
                   "EulerianGraphQ[Graph[{1,2,3,4},{1->2,2->3,3->1}]], "
                   "EulerianGraphQ[Graph[{1,2,3,4,5,6},{1<->2,2<->3,3<->1,4<->5,5<->6,6<->4}]], "
                   "EulerianGraphQ[Graph[{1},{}]], EulerianGraphQ[Graph[{},{}]], "
                   "EulerianGraphQ[Graph[{1,2},{}]], EulerianGraphQ[5]}",
                   "{True, False, True, True, False, True, False, True, False}", 0);
    assert_eval_eq("Head[EulerianGraphQ[Graph[{1,2,3},{1->2, 2<->3, 3->1}]]]", "EulerianGraphQ", 0);
}

static void test_eulerian_cycle(void) {
    assert_eval_eq("FindEulerianCycle[CycleGraph[4]]", "{{1 <-> 4, 4 <-> 3, 3 <-> 2, 2 <-> 1}}", 0);
    assert_eval_eq("FindEulerianCycle[" G4 "]", "{}", 0);
    assert_eval_eq("FindEulerianCycle[Graph[{1,2,3},{1->2,2->3,3->1}]]", "{{1 -> 2, 2 -> 3, 3 -> 1}}", 0);
    assert_eval_eq("FindEulerianCycle[Graph[{1,2,3,4,5},{1<->2,2<->3,3<->1,1<->4,4<->5,5<->1}]]",
                   "{{1 <-> 5, 5 <-> 4, 4 <-> 1, 1 <-> 3, 3 <-> 2, 2 <-> 1}}", 0);
    assert_eval_eq("FindEulerianCycle[Graph[{1,2,3,4,5,6},{1<->2,2<->3,3<->1,3<->4,4<->5,5<->3}]]",
                   "{{1 <-> 3, 3 <-> 5, 5 <-> 4, 4 <-> 3, 3 <-> 2, 2 <-> 1}}", 0);
    assert_eval_eq("FindEulerianCycle[Graph[{1,2,3,4,5},{1->2,2->3,3->1,1->4,4->5,5->1}]]",
                   "{{1 -> 2, 2 -> 3, 3 -> 1, 1 -> 4, 4 -> 5, 5 -> 1}}", 0);
    assert_eval_eq("FindEulerianCycle[Graph[{3,2,1},{1<->2,2<->3,3<->1}]]",
                   "{{3 <-> 1, 1 <-> 2, 2 <-> 3}}", 0);
    assert_eval_eq("{FindEulerianCycle[Graph[{1},{}]], FindEulerianCycle[Graph[{},{}]], "
                   "FindEulerianCycle[CycleGraph[5], 1]}",
                   "{{{}}, {}, {{1 <-> 5, 5 <-> 4, 4 <-> 3, 3 <-> 2, 2 <-> 1}}}", 0);
    /* a 4-regular circulant: every edge exactly once */
    assert_eval_eq("eu = Graph[Range[40], Join[Table[i <-> Mod[i,40]+1, {i,40}], "
                   "Table[i <-> Mod[i+1,40]+1, {i,40}]]]; c = First[FindEulerianCycle[eu]]; "
                   "{Length[c], Sort[Sort /@ (List @@@ c)] === Sort[Sort /@ (List @@@ EdgeList[eu])], "
                   "And @@ Table[c[[i, 2]] === c[[Mod[i, 80] + 1, 1]], {i, 80}]}",
                   "{80, True, True}", 0);
    assert_eval_eq("Head[FindEulerianCycle[CycleGraph[4], 2]]", "FindEulerianCycle", 0);
}

static void test_find_cycle(void) {
    assert_eval_eq("FindCycle[" G4 "]", "{}", 0);
    assert_eval_eq("FindCycle[CycleGraph[4]]", "{{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}}", 0);
    assert_eval_eq("FindCycle[Graph[{1<->2,2<->3,3<->4,4<->2}]]", "{{2 <-> 3, 3 <-> 4, 4 <-> 2}}", 0);
    assert_eval_eq("FindCycle[Graph[{1<->2,2<->3,3<->1,3<->4,4<->5,5<->3}]]",
                   "{{3 <-> 4, 4 <-> 5, 5 <-> 3}}", 0);
    assert_eval_eq("FindCycle[Graph[{1->2,2->3,3->4,4->2,3->1}]]", "{{1 -> 2, 2 -> 3, 3 -> 1}}", 0);
    assert_eval_eq("FindCycle[Graph[{1,2},{1->2,2->1}]]", "{{1 -> 2, 2 -> 1}}", 0);
    /* length-bounded and enumerating forms (Mathilda's own order) */
    assert_eval_eq("FindCycle[CompleteGraph[4], 3]", "{{1 <-> 2, 2 <-> 3, 3 <-> 1}}", 0);
    assert_eval_eq("Length[FindCycle[CompleteGraph[4], {3}, All]]", "4", 0);
    assert_eval_eq("Length[FindCycle[CompleteGraph[4], Infinity, All]]", "7", 0);
    assert_eval_eq("Length[FindCycle[CompleteGraph[6], Infinity, All]]", "197", 0);
    assert_eval_eq("FindCycle[CompleteGraph[4], {4}]", "{{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}}", 0);
    assert_eval_eq("FindCycle[Graph[{1->2,2->3,3->1}], {2}]", "{}", 0);
    assert_eval_eq("FindCycle[Graph[{1->2,2->1,2->3,3->1}], Infinity, All]",
                   "{{1 -> 2, 2 -> 1}, {1 -> 2, 2 -> 3, 3 -> 1}}", 0);
    assert_eval_eq("FindCycle[{CompleteGraph[4], 2}, 3, All]",
                   "{{2 <-> 1, 1 <-> 3, 3 <-> 2}, {2 <-> 1, 1 <-> 4, 4 <-> 2}, {2 <-> 3, 3 <-> 4, 4 <-> 2}}", 0);
    assert_eval_eq("FindCycle[{Graph[{1->2,2->3,3->1,3->4,4->2}], 4}]", "{{4 -> 2, 2 -> 3, 3 -> 4}}", 0);
    assert_eval_eq("FindCycle[{Graph[{1<->2,2<->3,3<->4,4<->5,5<->2}], 1}]", "{}", 0);
    assert_eval_eq("TimeConstrained[FindCycle[CompleteGraph[30], {30}, All], 0.3, tout]", "tout", 0);
    assert_eval_eq("Head[FindCycle[Graph[{1,2,3},{1->2, 2<->3, 3->1}]]]", "FindCycle", 0);
}

static void test_find_path(void) {
    assert_eval_eq("FindPath[" G4 ", 1, 4]", "{{1, 2, 3, 4}}", 0);
    assert_eval_eq("FindPath[CompleteGraph[4], 1, 4]", "{{1, 2, 3, 4}}", 0);
    assert_eval_eq("FindPath[" D3 ", 3, 1]", "{}", 0);
    assert_eval_eq("FindPath[CompleteGraph[4], 1, 4, 2, All]", "{{1, 4}, {1, 3, 4}, {1, 2, 4}}", 0);
    assert_eval_eq("FindPath[CompleteGraph[4], 1, 4, {2}, All]", "{{1, 3, 4}, {1, 2, 4}}", 0);
    assert_eval_eq("FindPath[CompleteGraph[4], 1, 4, 3, 2]", "{{1, 2, 4}, {1, 2, 3, 4}}", 0);
    assert_eval_eq("FindPath[Graph[{1,2,3,4,5},{1<->2,2<->3,1<->4,4<->5,5<->3}], 1, 3, Infinity, All]",
                   "{{1, 2, 3}, {1, 4, 5, 3}}", 0);
    assert_eval_eq("FindPath[Graph[{1,2,3,4,5},{1<->4,1<->2,2<->3,4<->5,5<->3}], 1, 3]",
                   "{{1, 4, 5, 3}}", 0);
    assert_eval_eq("FindPath[CompleteGraph[4], 1, 1]", "{}", 0);
    assert_eval_eq("Length[FindPath[CompleteGraph[6], 1, 6, Infinity, All]]", "65", 0);
    assert_eval_eq("Head[FindPath[CompleteGraph[4], 1, 9]]", "FindPath", 0);
}

static void test_memo_and_results(void) {
    /* results are valid, memoized graphs usable by every accessor */
    assert_eval_eq("h = VertexDelete[PathGraph[Range[6]], 3]; {GraphQ[h], ConnectedGraphQ[h], "
                   "EdgeQ[h, 4<->5], VertexQ[h, 3]}", "{True, False, True, False}", 0);
    assert_eval_eq("h = GraphUnion[CycleGraph[5], Graph[{1<->3}]]; {EdgeCount[h], EdgeQ[h, 3<->1]}",
                   "{6, True}", 0);
    assert_eval_eq("InputForm[EdgeAdd[Graph[{1,2},{1<->2}], 2<->3]]",
                   "Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}]", 0);
    /* non-graph arguments stay unevaluated */
    assert_eval_eq("{Head[VertexAdd[5, 1]], Head[GraphUnion[5, 6]], Head[FindPath[x, 1, 2]], "
                   "Head[Subgraph[x, {1}]]}", "{VertexAdd, GraphUnion, FindPath, Subgraph}", 0);
}

int main(void) {
    symtab_init();
    core_init();
    setup();

    TEST(test_vertex_add_delete);
    TEST(test_edge_add_delete);
    TEST(test_subgraph_neighborhood);
    TEST(test_replace_rules_index);
    TEST(test_transforms);
    TEST(test_set_operations);
    TEST(test_predicates);
    TEST(test_eulerian_cycle);
    TEST(test_find_cycle);
    TEST(test_find_path);
    TEST(test_memo_and_results);

    printf("All graph ops tests passed!\n");
    return 0;
}
