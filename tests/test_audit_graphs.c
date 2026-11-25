/* Builtin-audit regressions: graphs (review/issues/<Head>.md; tests/test_audit.h). */
#include "test_audit.h"

/* CompleteGraph#1 (critical): CompleteGraph[10^7] crashed (unchecked calloc
 * of n(n-1)/2 edge slots).  The documented resource limit (10^8 vertices or
 * 5x10^7 edges) leaves such a family unevaluated. */
static void test_complete_graph_resource_limit(void) {
    audit_chk("CompleteGraph#1", "CompleteGraph[10^7]", "CompleteGraph[10000000]");
    audit_chk("CompleteGraph#1", "CompleteGraph[2^62]", "CompleteGraph[4611686018427387904]");
    audit_chk("CompleteGraph#1", "EdgeCount[CompleteGraph[5]]", "10");
    audit_chk("CompleteGraph#1", "EdgeCount[CompleteGraph[1]]", "0");
    audit_chk("CompleteGraph#1", "VertexCount[CompleteGraph[0]]", "0");
}

/* StarGraph#1 (high): StarGraph[10^15] crashed (unchecked calloc). */
static void test_star_graph_resource_limit(void) {
    audit_chk("StarGraph#1", "StarGraph[10^15]", "StarGraph[1000000000000000]");
    audit_chk("StarGraph#1", "EdgeCount[StarGraph[6]]", "5");
    audit_chk("StarGraph#1", "CycleGraph[10^15]", "CycleGraph[1000000000000000]");
    audit_chk("StarGraph#1", "PathGraph[10^15]", "PathGraph[1000000000000000]");
    audit_chk("StarGraph#1", "RandomGraph[{2^40, 1}]", "RandomGraph[{1099511627776, 1}]");
    audit_chk("StarGraph#1", "EdgeCount[CycleGraph[4]]", "4");
}

/* FindSpanningTree#1 (critical) / #2: "Weighted directed FindSpanningTree
 * returns a subgraph containing directed cycles (not a branching)", and the
 * result depended on earlier evaluations.  Chu-Liu/Edmonds recorded every
 * edge internal to a contracted super-vertex as a one-vertex cycle and wrote
 * past the N-entry cycle arrays.  The issue's 12-vertex graph has minimum
 * arborescence weight 28; random small digraphs are checked against an
 * exhaustive enumeration of every root's parent choices. */
static void test_find_spanning_tree_directed(void) {
    audit_run("fstWs = {11,2,-5,19,1,12,12,2,7,11,6,13,6,9,3,16,12,14,"
              "18,-5,7,20,18,11,20,-1,11,19,12,1,8,-4,10,6,13,12,1,11,8,10}");
    audit_run("fstG = Graph[Range[12],{1->7,1->12,2->3,2->5,2->6,2->8,3->6,3->10,3->12,4->2,"
              "5->4,5->10,5->12,6->1,6->4,6->8,7->1,7->9,7->11,7->12,8->1,8->4,8->5,8->12,"
              "9->4,9->7,9->11,10->1,10->2,10->8,10->12,11->4,11->5,11->7,11->9,12->1,12->6,"
              "12->7,12->8,12->9},EdgeWeight->fstWs]");
    /* tree weight, read back through the input edge list and weights */
    audit_run("fstW[g_, ws_, t_] := Module[{el = EdgeList[g]},"
              " Total[ws[[Flatten[Position[el, #, {1}, 1] & /@ EdgeList[t]]]]]]");
    audit_run("fstArbQ[n_, t_] := EdgeCount[t] == n - 1 && AcyclicGraphQ[t] &&"
              " Max[Counts[EdgeList[t][[All, 2]]]] == 1");
    audit_run("fstT = FindSpanningTree[fstG]");
    audit_true("FindSpanningTree#1", "fstArbQ[12, fstT]");
    audit_chk("FindSpanningTree#1", "fstW[fstG, fstWs, fstT]", "28");
    audit_chk("FindSpanningTree#1", "FindSpanningTree[fstG] === fstT", "True");
    audit_run("fstWs4 = {7,1,2,9,6,6,6,8,8,2,2,5}; fstG4 = Graph[Range[4],{1->2,1->3,1->4,"
              "2->1,2->3,2->4,3->1,3->2,3->4,4->1,4->2,4->3},EdgeWeight->fstWs4]");
    audit_true("FindSpanningTree#2", "fstArbQ[4, FindSpanningTree[fstG4]]");
    audit_chk("FindSpanningTree#2", "fstW[fstG4, fstWs4, FindSpanningTree[fstG4]]", "5");
    audit_chk("FindSpanningTree#2", "FindSpanningTree[fstG4] === FindSpanningTree[fstG4]", "True");
    /* brute force: every root r, every choice of one incoming edge per other
     * vertex; a choice is an arborescence when following parents n times
     * from every vertex lands on r. */
    audit_run("fstBrute[n_, edges_, ws_] := Module[{best = Infinity, inc, others, p, w},"
              " Do[inc = Table[Select[Range[Length[edges]], edges[[#, 2]] == v &], {v, n}];"
              "  others = Delete[Range[n], r];"
              "  If[AllTrue[inc[[others]], Length[#] > 0 &],"
              "   Do[p = Range[n]; w = 0;"
              "    Do[p[[others[[k]]]] = edges[[ch[[k]], 1]]; w += ws[[ch[[k]]]], {k, Length[others]}];"
              "    If[Nest[#[[#]] &, p, n] === ConstantArray[r, n] && w < best, best = w],"
              "    {ch, Tuples[inc[[others]]]}]], {r, n}]; best]");
    audit_run("SeedRandom[1234]; fstBad = 0; fstCnt = 0;"
              " Do[n = RandomInteger[{3, 5}];"
              "  edges = Select[Tuples[Range[n], 2], #[[1]] != #[[2]] && RandomReal[] < 0.7 &];"
              "  ws = RandomInteger[{-6, 12}, Length[edges]];"
              "  b = fstBrute[n, edges, ws];"
              "  If[b =!= Infinity, fstCnt++;"
              "   g = Graph[Range[n], Rule @@@ edges, EdgeWeight -> ws]; t = FindSpanningTree[g];"
              "   If[!(fstArbQ[n, t] && fstW[g, ws, t] == b), fstBad++]], {40}]");
    audit_true("FindSpanningTree#1", "fstCnt >= 30");
    audit_chk("FindSpanningTree#1", "fstBad", "0");
}

/* CirculantGraph#1 (low): the 5x10^7 edge cap counted raw edges before
 * deduplication, so a repeated or half-turn jump was refused far below the
 * limit.  Jumps are now canonicalised (mod n, j ~ n - j, distinct, the
 * half-turn emitted once), so every emitted edge is distinct.  The large case
 * (CirculantGraph[2*10^7, {1, 1, 1}] -> 2*10^7 edges) takes ~5 GB, so it is
 * checked by hand; the canonicalisation is pinned here on small n. */
static void test_circulant_graph_dedup_before_cap(void) {
    audit_chk("CirculantGraph#1", "EdgeCount[CirculantGraph[10, {1, 1, 1}]]", "10");
    audit_chk("CirculantGraph#1", "EdgeCount[CirculantGraph[10, {1, 9, -1, 11}]]", "10");
    audit_chk("CirculantGraph#1", "EdgeCount[CirculantGraph[10, 5]]", "5");
    audit_chk("CirculantGraph#1", "EdgeCount[CirculantGraph[10, {1, 5}]]", "15");
    audit_chk("CirculantGraph#1", "EdgeCount[CirculantGraph[10, {0, 10}]]", "0");
    audit_true("CirculantGraph#1",
               "EdgeList[CirculantGraph[6, {1, 3, 1}]] === EdgeList[CirculantGraph[6, {1, 3}]]");
    audit_chk("CirculantGraph#1", "EdgeList[CirculantGraph[6, 2]]",
              "{1 <-> 3, 1 <-> 5, 2 <-> 4, 2 <-> 6, 3 <-> 5, 4 <-> 6}");
    audit_chk("CirculantGraph#1", "EdgeList[CirculantGraph[4, 2]]", "{1 <-> 3, 2 <-> 4}");
}

/* FindMaximumFlow#1 (high): "Small capacities next to a very large one are
 * rounded to 0, giving max flow 0.0".  On a path the maximum flow is the
 * bottleneck capacity. */
static void test_max_flow_wide_range(void) {
    audit_chk("FindMaximumFlow#1", "FindMaximumFlow[Graph[{1 -> 2, 2 -> 3}], 1, 3, EdgeCapacity -> {1.*^20, 1.}]", "1.0");
    audit_chk("FindMaximumFlow#1", "FindMaximumFlow[Graph[{1 -> 2, 2 -> 3}], 1, 3, EdgeCapacity -> {1.*^20, 3.}]", "3.0");
    audit_chk("FindMaximumFlow#1", "FindMaximumFlow[Graph[{1 -> 2, 2 -> 3}], 1, 3, EdgeCapacity -> {2^63, 5}]", "5");
    audit_chk("FindMaximumFlow#1", "FindMaximumFlow[Graph[{1 -> 2, 2 -> 3}], 1, 3, VertexCapacity -> {1.*^20, 1., 5.}]", "1.0");
    audit_chk("FindMaximumFlow#1", "FindMaximumFlow[Graph[{1 -> 2, 2 -> 3}], 1, 3, EdgeCapacity -> {1.*^300, 1.*^-300}]", "1.0*^-300");
    audit_chk("FindMaximumFlow#1", "FindMaximumFlow[Graph[{1 -> 2, 2 -> 3}], 1, 3, EdgeCapacity -> {Infinity, 1.*^-300}]", "1.0*^-300");
    audit_chk("FindMaximumFlow#1", "FindMaximumFlow[Graph[{1 -> 2, 2 -> 3}], 1, 3, \"FlowMatrix\", EdgeCapacity -> {1.*^20, 1.}]",
              "{{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, 0.0}}");
    /* an integer network whose capacity sum overflows int64 is exact while
     * the flow itself fits (min cut {7, 10^18}) */
    audit_chk("FindMaximumFlow#1", "FindMaximumFlow[Graph[{1 -> 2, 2 -> 3, 1 -> 3}], 1, 3, EdgeCapacity -> {10^30, 7, 10^18}]",
              "1000000000000000007");
    /* against the same network with the small capacities scaled to integers
     * (exact) and the large ones lowered (they are not on the minimum cut) */
    audit_run("SeedRandom[7]; mfg = RandomGraph[{30, 120}]; "
              "mfd = Graph[VertexList[mfg], DirectedEdge @@@ EdgeList[mfg]]; "
              "mfc = RandomChoice[{1., 1.*^18, 3.5, 1.*^-5}, EdgeCount[mfd]]");
    audit_true("FindMaximumFlow#1",
               "Abs[FindMaximumFlow[mfd, 1, 30, EdgeCapacity -> mfc] - "
               "FindMaximumFlow[mfd, 1, 30, EdgeCapacity -> Round[10^5 (mfc /. 1.*^18 -> 10.^6)]]/10^5] < 10^-12");
    /* sanity */
    audit_chk("FindMaximumFlow#1", "FindMaximumFlow[Graph[{1 -> 2, 2 -> 3}], 1, 3, EdgeCapacity -> {3, 5}]", "3");
    audit_chk("FindMaximumFlow#1", "FindMaximumFlow[Graph[{1 -> 2, 2 -> 3}], 1, 3, EdgeCapacity -> {0.5, 2.25}]", "0.5");
    audit_chk("FindMaximumFlow#1", "FindMaximumFlow[Graph[{1 -> 2, 2 -> 3}], 1, 3, EdgeCapacity -> {Infinity, Infinity}]", "Infinity");
}

/* FindMinimumCut#1 (high): "FindMinimumCut returns cut value 0.0 when real
 * weights span a wide range"; #2 is the same rounding for a tiny weight. */
static void test_min_cut_wide_range(void) {
    audit_chk("FindMinimumCut#1", "FindMinimumCut[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}, EdgeWeight -> {1.*^20, 1.}]]",
              "{1.0, {{3}, {1, 2}}}");
    audit_chk("FindMinimumCut#1", "FindMinimumCut[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}, EdgeWeight -> {1.*^-20, 3.}]]",
              "{1.0*^-20, {{2, 3}, {1}}}");
    audit_chk("FindMinimumCut#1", "FindMinimumCut[Graph[{1, 2, 3}, {1 -> 2, 2 -> 3, 3 -> 1}, EdgeWeight -> {1.*^20, 1., 1.*^20}]][[1]]",
              "1.0");
    audit_chk("FindMinimumCut#1", "EdgeConnectivity[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}, EdgeWeight -> {1.*^20, 1.}]]", "1.0");
    audit_chk("FindMinimumCut#1", "FindEdgeCut[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}, EdgeWeight -> {1.*^20, 1.}], 1, 3]", "{2 <-> 3}");
    audit_run("SeedRandom[7]; mcg = RandomGraph[{30, 120}]; mcw = RandomChoice[{1., 1.*^18, 3.5, 1.*^-5}, EdgeCount[mcg]]");
    audit_true("FindMinimumCut#1",
               "Abs[FindMinimumCut[Graph[VertexList[mcg], EdgeList[mcg], EdgeWeight -> mcw]][[1]] - "
               "FindMinimumCut[Graph[VertexList[mcg], EdgeList[mcg], EdgeWeight -> Round[10^5 (mcw /. 1.*^18 -> 10.^6)]]][[1]]/10^5] < 10^-12");
    /* sanity */
    audit_chk("FindMinimumCut#1", "FindMinimumCut[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}, EdgeWeight -> {2, 5}]]",
              "{2, {{2, 3}, {1}}}");
}

/* FindMinimumCut#3-4, FindMaximumFlow#3, EdgeConnectivity#2: exact Integer
 * capacities past 2^61 give exact answers (they saturated at 2^61 or came
 * back Real). */
static void test_s14_exact_big_capacities(void) {
    audit_true("FindMinimumCut#3",
               "FindMinimumCut[Graph[{1, 2}, {1 <-> 2}, EdgeWeight -> {10^20}]] === {10^20, {{2}, {1}}}");
    audit_true("FindMinimumCut#4",
               "FindMinimumCut[Graph[{1, 2}, {1 <-> 2}, EdgeWeight -> {2^62}]] === {2^62, {{2}, {1}}}");
    audit_true("FindMinimumCut#4",
               "First[FindMinimumCut[Graph[{1, 2}, {1 <-> 2}, EdgeWeight -> {2^61 + 3}]]] === 2^61 + 3");
    audit_true("FindMaximumFlow#3",
               "FindMaximumFlow[Graph[{1 -> 2, 2 -> 3}], 1, 3, EdgeCapacity -> {2^61, 7}] === 7");
    audit_true("FindMaximumFlow#3", "FindMaximumFlow[Graph[{1 -> 2, 2 -> 3, 1 -> 3}], 1, 3, "
                                    "EdgeCapacity -> {2^62, 2^62, 2^62}] === 2^63");
    audit_true("EdgeConnectivity#2", "EdgeConnectivity[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 1}, "
                                     "EdgeWeight -> {2^61, 2^61, 2^61}]] === 2^62");
    /* a scaled weighted graph: every cut and flow scales exactly */
    audit_run("SeedRandom[3]; auditG = RandomGraph[{40, 200}]; "
              "auditW = RandomInteger[{1, 50}, 200]");
    audit_true("FindMinimumCut#3",
               "First[FindMinimumCut[Graph[EdgeList[auditG], EdgeWeight -> 3^50 auditW]]] === "
               "3^50 First[FindMinimumCut[Graph[EdgeList[auditG], EdgeWeight -> auditW]]]");
    audit_true("FindMaximumFlow#3",
               "FindMaximumFlow[Graph[EdgeList[auditG]], 1, 2, EdgeCapacity -> 3^50 auditW] === "
               "3^50 FindMaximumFlow[Graph[EdgeList[auditG]], 1, 2, EdgeCapacity -> auditW]");
    audit_run("ClearAll[auditG, auditW]");
}

/* Wave 4: generator options, duplicate vertices, pattern / edge-list forms,
 * integer truncation, weights and capacities. */
static void test_wave4_graphs(void) {
    /* generators take DirectedEdges / edge options (graph_generator_with_options) */
    audit_chk("CycleGraph#1", "EdgeList[CycleGraph[5, DirectedEdges -> True]]",
              "{1 -> 2, 2 -> 3, 3 -> 4, 4 -> 5, 5 -> 1}");
    audit_chk("PathGraph#2", "EdgeList[PathGraph[{1, 2, 3}, DirectedEdges -> True]]",
              "{1 -> 2, 2 -> 3}");
    audit_chk("StarGraph#2", "EdgeList[StarGraph[4, DirectedEdges -> True]]",
              "{1 -> 2, 1 -> 3, 1 -> 4}");
    audit_chk("AdjacencyGraph#3", "EdgeList[AdjacencyGraph[{{0, 1}, {1, 0}}, DirectedEdges -> True]]",
              "{1 -> 2, 2 -> 1}");
    audit_chk("AdjacencyGraph#2", "EdgeList[AdjacencyGraph[{a, b, c}, {{0,1,0},{1,0,1},{0,1,0}}]]",
              "{a <-> b, b <-> c}");
    /* entries are read as full integers: only 0 and 1 are edges */
    audit_true("AdjacencyGraph#1",
               "Head[AdjacencyGraph[{{0, 4294967297}, {4294967297, 0}}]] === AdjacencyGraph");
    audit_true("AdjacencyGraph#1",
               "Head[AdjacencyGraph[{{0, 4294967296}, {4294967296, 0}}]] === AdjacencyGraph");
    /* a repeated vertex is one vertex */
    audit_chk("Graph#1", "g = Graph[{1, 1, 2}, {1 -> 2}]; {GraphQ[g], VertexList[g], AdjacencyMatrix[g]}",
              "{True, {1, 2}, {{0, 1}, {0, 0}}}");
    audit_chk("VertexCount#2", "VertexCount[Graph[{1, 1, 2}, {1 -> 2}]]", "2");
    audit_chk("PathGraph#1", "{VertexCount[PathGraph[{a, b, c, a}]], VertexDegree[PathGraph[{a, b, c, a}]]}",
              "{3, {2, 2, 2}}");
    audit_run("ClearAll[g]");
    /* mixed graph: parallel directed + undirected edges are counted */
    audit_chk("AdjacencyMatrix#1", "AdjacencyMatrix[Graph[{1 -> 2, 2 -> 1, 1 <-> 2}]]", "{{0, 2}, {2, 0}}");
    /* pattern forms */
    audit_chk("VertexCount#1", "{VertexCount[Graph[{1, 2, 3}, {1 -> 2}], 1 | 2], "
              "VertexCount[Graph[{a, b, c}, {}], _Symbol]}", "{2, 3}");
    audit_chk("EdgeCount#1", "EdgeCount[Graph[{1 -> 2, 2 -> 3, 1 -> 3}], 1 -> _]", "2");
    audit_chk("EdgeList#1", "EdgeList[Graph[{1 -> 2, 2 -> 3}], DirectedEdge[1, _]]", "{1 -> 2}");
    audit_chk("VertexList#1", "VertexList[Graph[{1 -> 2, 2 -> 3, 3 -> 1}], _?OddQ]", "{1, 3}");
    /* a list of edges stands for its graph */
    audit_chk("VertexList#2", "VertexList[{1 -> 2, 2 -> 3}]", "{1, 2, 3}");
    audit_chk("ConnectedComponents#1", "ConnectedComponents[{1 <-> 2, 3 <-> 4, 4 <-> 5}]",
              "{{3, 4, 5}, {1, 2}}");
    audit_chk("ConnectedComponents#1", "ConnectedComponents[{1 -> 2, 2 -> 1, 3 -> 2}]", "{{1, 2}, {3}}");
    /* vertices within distance d */
    audit_chk("AdjacencyList#1", "AdjacencyList[Graph[{1, 2, 3, 4}, {1 <-> 2, 2 <-> 3, 3 <-> 4}], 1, 2]",
              "{2, 3}");
    /* s-t vertex connectivity */
    audit_chk("VertexConnectivity#1", "{VertexConnectivity[CycleGraph[5], 1, 3], "
              "VertexConnectivity[PathGraph[{1, 2, 3, 4}], 1, 4]}", "{2, 1}");
    /* size bounds saturate instead of truncating to int */
    audit_chk("FindIndependentVertexSet#1", "FindIndependentVertexSet[CycleGraph[6], 4294967296]", "{{1, 3, 5}}");
    audit_chk("FindIndependentVertexSet#1", "FindIndependentVertexSet[CycleGraph[6], 4294967298]", "{{1, 3, 5}}");
    /* large k: a star, no overflow */
    audit_chk("KaryTree#1", "EdgeList[KaryTree[3, 2^62]]", "{1 <-> 2, 1 <-> 3}");
    audit_chk("KaryTree#1", "EdgeList[KaryTree[5, 2^61]]",
              "{1 <-> 2, 1 <-> 3, 1 <-> 4, 1 <-> 5}");
    /* numeric-constant weights are weights */
    audit_true("GraphDistance#1", "Abs[GraphDistance[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3, 1 <-> 3}, "
               "EdgeWeight -> {Pi, 1, 10}], 1, 3] - N[Pi + 1]] < 1*^-12");
    /* EdgeAdd keeps EdgeCapacity */
    audit_chk("EdgeAdd#1", "FindMaximumFlow[EdgeAdd[Graph[{1 -> 2}, EdgeCapacity -> {7}], 2 -> 3], 1, 2]", "7");
    /* periodic chain at a = 1 */
    audit_true("PageRankCentrality#1",
               "Max[Abs[PageRankCentrality[PathGraph[3], 1] - {0.25, 0.5, 0.25}]] < 1*^-9");
    /* the n shortest paths, shortest first */
    audit_true("FindPath#1", "With[{p = FindPath[CompleteGraph[4], 1, 4, Infinity, 2]}, "
               "First[p] === {1, 4} && Length[p] == 2 && Length[p[[2]]] == 3]");
    /* past the 8192-vertex bitset engine */
    audit_true("FindClique#1", "With[{c = FindClique[CycleGraph[10000], 2]}, "
               "Head[c] === List && Length[c] == 1 && Length[First[c]] == 2]");
    audit_true("FindKClique#1", "With[{r = FindKClique[PathGraph[Range[10000]], 2]}, "
               "Head[r] === List && Length[First[r]] == 3]");
    /* an unbalanced bipartite graph has no Hamiltonian cycle */
    audit_chk("HamiltonianGraphQ#1", "HamiltonianGraphQ[CompleteGraph[{8, 9}]]", "False");
    audit_chk("HamiltonianGraphQ#2", "TimeConstrained[HamiltonianGraphQ[GridGraph[{11, 11}]], 10]", "False");
}

int main(void) {
    audit_init();
    TEST(test_wave4_graphs);
    TEST(test_complete_graph_resource_limit);
    TEST(test_star_graph_resource_limit);
    TEST(test_find_spanning_tree_directed);
    TEST(test_circulant_graph_dedup_before_cap);
    TEST(test_max_flow_wide_range);
    TEST(test_min_cut_wide_range);
    TEST(test_s14_exact_big_capacities);
    return audit_finish("audit_graphs");
}
