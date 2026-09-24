/* test_hypergraph.c - the native hypergraph subsystem (src/graph/hyp_*.c).
 *
 * Construction/validation/printing and the InputForm round-trip, the shared
 * accessors (VertexList ... IncidenceMatrix) on Hypergraph arguments, every
 * structural operation on small hand-checked cases, the memo's soundness under
 * in-place mutation, RandomHypergraph's shape and seeding, and transversals
 * cross-checked against brute force on random instances.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "print.h"
#include "graph.h"
#include "graph_hyper.h"
#include "test_utils.h"
#include <stdlib.h>

#define H1 "Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]"
#define H2 "Hypergraph[{{1,2,3},{2,3,4},{3,4,5},{1,5}}]"

static void test_construction(void) {
    assert_eval_eq(H1, "Hypergraph[<7 vertices, 4 hyperedges>]", 0);
    assert_eval_eq("Hypergraph[{1},{{1}}]", "Hypergraph[<1 vertex, 1 hyperedge>]", 0);
    assert_eval_eq("Hypergraph[{},{}]", "Hypergraph[<0 vertices, 0 hyperedges>]", 0);
    assert_eval_eq("InputForm[" H1 "]",
        "Hypergraph[{1, 2, 3, 4, 5, 6, 7}, {{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]", 0);
    assert_eval_eq(H1, "Hypergraph[List[1, 2, 3, 4, 5, 6, 7], "
        "List[List[1, 2, 3], List[3, 4], List[4, 5, 6], List[7]]]", 1);
    /* Round trip through the parser. */
    assert_eval_eq("ToExpression[ToString[InputForm[" H1 "]]] === " H1, "True", 0);
    /* Explicit vertices: order kept, duplicates dropped, isolated allowed. */
    assert_eval_eq("VertexList[Hypergraph[{3,1,2,1,9},{{1,2}}]]", "{3, 1, 2, 9}", 0);
    /* Repeated hyperedges, repeated vertices inside one, empty hyperedges. */
    assert_eval_eq("Hypergraph[{{1,1,2},{1,1,2},{}}]", "Hypergraph[<2 vertices, 3 hyperedges>]", 0);
    /* Integer fast path vs. the hash path: a Real is never SameQ an Integer
     * vertex; a huge value range and mixed vertex types take the hash index. */
    assert_eval_eq("HypergraphQ[Hypergraph[{1,2},{{1.,2}}]]", "False", 0);
    assert_eval_eq("VertexList[Hypergraph[{{10^12, 1}, {1, -5}}]]", "{1000000000000, 1, -5}", 0);
    assert_eval_eq("VertexDegree[Hypergraph[{{10^12, 1}, {1, -5}}], 10^12]", "1", 0);
    assert_eval_eq("VertexList[Hypergraph[{{a, 1}, {1, b}, {-3}}]]", "{a, 1, b, -3}", 0);
    assert_eval_eq("VertexDegree[Hypergraph[{-2, 0, 7}, {{7, -2}, {0, 7}}]]", "{1, 1, 2}", 0);
    /* Graph conversion. */
    assert_eval_eq("EdgeList[Hypergraph[CycleGraph[3]]]", "{{1, 2}, {2, 3}, {3, 1}}", 0);
    /* Malformed: unknown vertex, non-List hyperedge, wrong arity. */
    assert_eval_eq("HypergraphQ[Hypergraph[{1,2},{{1,5}}]]", "False", 0);
    assert_eval_eq("HypergraphQ[Hypergraph[{{1,2},3}]]", "False", 0);
    assert_eval_eq("HypergraphQ[Hypergraph[{1},{{1}},x]]", "False", 0);
    assert_eval_eq("{HypergraphQ[" H1 "], HypergraphQ[5], HypergraphQ[CycleGraph[3]]}",
                   "{True, False, False}", 0);
    /* Graph heads keep rejecting hypergraphs and vice versa. */
    assert_eval_eq("GraphQ[" H1 "]", "False", 0);
}

static void test_accessors(void) {
    assert_eval_eq("VertexList[" H1 "]", "{1, 2, 3, 4, 5, 6, 7}", 0);
    assert_eval_eq("EdgeList[" H1 "]", "{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}", 0);
    assert_eval_eq("{VertexCount[" H1 "], EdgeCount[" H1 "]}", "{7, 4}", 0);
    assert_eval_eq("VertexDegree[" H1 "]", "{1, 1, 2, 2, 1, 1, 1}", 0);
    assert_eval_eq("VertexDegree[" H1 ", 4]", "2", 0);
    assert_eval_eq("VertexDegree[" H1 ", 99]", "VertexDegree[Hypergraph[<7 vertices, 4 hyperedges>], 99]", 0);
    /* A vertex repeated in one hyperedge counts once; repeated hyperedges count. */
    assert_eval_eq("VertexDegree[Hypergraph[{{1,1,2},{1,1,2}}]]", "{2, 2}", 0);
    assert_eval_eq("IncidenceMatrix[Hypergraph[{{1,2},{2,3,3}}]]", "{{1, 0}, {1, 1}, {0, 2}}", 0);
    assert_eval_eq("IncidenceMatrix[" H1 "]",
        "{{1, 0, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {0, 1, 1, 0}, {0, 0, 1, 0}, "
        "{0, 0, 1, 0}, {0, 0, 0, 1}}", 0);
    /* Directed degree heads are Graph-only. */
    assert_eval_eq("Head[VertexInDegree[" H1 "]]", "VertexInDegree", 0);
    /* Graph behaviour unchanged. */
    assert_eval_eq("VertexDegree[Graph[{1->2,2->3}]]", "{1, 2, 1}", 0);
    assert_eval_eq("IncidenceMatrix[Graph[{1->2}]]", "{{-1}, {1}}", 0);
}

static void test_sizes(void) {
    assert_eval_eq("HyperedgeSizes[" H1 "]", "{3, 2, 3, 1}", 0);
    assert_eval_eq("{HypergraphRank[" H1 "], HypergraphCorank[" H1 "]}", "{3, 1}", 0);
    assert_eval_eq("{HypergraphRank[Hypergraph[{1},{}]], HypergraphCorank[Hypergraph[{1},{}]]}",
                   "{0, 0}", 0);
    assert_eval_eq("{UniformHypergraphQ[" H1 "], UniformHypergraphQ[" H2 "], "
                   "UniformHypergraphQ[Hypergraph[{{1,2,3},{2,3,4}}], 3], "
                   "UniformHypergraphQ[Hypergraph[{{1,2,3}}], 2], UniformHypergraphQ[7]}",
                   "{False, False, True, False, False}", 0);
}

static void test_dual_and_expansions(void) {
    assert_eval_eq("InputForm[HypergraphDual[" H1 "]]",
        "Hypergraph[{1, 2, 3, 4}, {{1}, {1}, {1, 2}, {2, 3}, {3}, {3}, {4}}]", 0);
    /* Dual of the dual recovers the incidence structure. */
    assert_eval_eq("EdgeList[HypergraphDual[HypergraphDual[" H2 "]]]",
                   "{{1, 2, 3}, {2, 3, 4}, {3, 4, 5}, {1, 5}}", 0);
    assert_eval_eq("InputForm[HypergraphCliqueExpansion[" H1 "]]",
        "Graph[{1, 2, 3, 4, 5, 6, 7}, {1 <-> 2, 1 <-> 3, 2 <-> 3, 3 <-> 4, "
        "4 <-> 5, 4 <-> 6, 5 <-> 6}]", 0);
    assert_eval_eq("InputForm[HypergraphStarExpansion[Hypergraph[{{1,2},{2,3}}]]]",
        "Graph[{1, 2, 3, Hyperedge[1], Hyperedge[2]}, {1 <-> Hyperedge[1], "
        "2 <-> Hyperedge[1], 2 <-> Hyperedge[2], 3 <-> Hyperedge[2]}]", 0);
    assert_eval_eq("BipartiteGraphQ[HypergraphStarExpansion[" H2 "]]", "True", 0);
    /* A vertex that is itself a Hyperedge[j] node: unevaluated. */
    assert_eval_eq("Head[HypergraphStarExpansion[Hypergraph[{{Hyperedge[1], 2}}]]]",
                   "HypergraphStarExpansion", 0);
    /* FR HypergraphToGraph: ordered, loops dropped, parallels merged. */
    assert_eval_eq("InputForm[HypergraphToGraph[{{1,2,3},{3,4}}]]",
        "Graph[{1, 2, 3, 4}, {1 -> 2, 1 -> 3, 2 -> 3, 3 -> 4}]", 0);
    assert_eval_eq("EdgeList[HypergraphToGraph[Hypergraph[{{1,1,2},{1,2},{2,1}}]]]",
                   "{1 -> 2, 2 -> 1}", 0);
}

static void test_line_graph(void) {
    assert_eval_eq("EdgeList[HypergraphLineGraph[" H1 "]]", "{1 <-> 2, 2 <-> 3}", 0);
    assert_eval_eq("EdgeList[HypergraphLineGraph[" H2 "]]",
                   "{1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 3, 3 <-> 4}", 0);
    assert_eval_eq("EdgeList[HypergraphLineGraph[" H2 ", 2]]", "{1 <-> 2, 2 <-> 3}", 0);
    assert_eval_eq("EdgeCount[HypergraphLineGraph[" H2 ", 3]]", "0", 0);
    assert_eval_eq("VertexList[HypergraphLineGraph[" H2 ", 3]]", "{1, 2, 3, 4}", 0);
    assert_eval_eq("Head[HypergraphLineGraph[" H2 ", 0]]", "HypergraphLineGraph", 0);
}

static void test_components_and_distances(void) {
    assert_eval_eq("HypergraphConnectedComponents[" H1 "]", "{{1, 2, 3, 4, 5, 6}, {7}}", 0);
    assert_eval_eq("HypergraphConnectedComponents[Hypergraph[{5,1,2,3},{{2,3}}]]",
                   "{{5}, {1}, {2, 3}}", 0);
    assert_eval_eq("{ConnectedHypergraphQ[" H1 "], ConnectedHypergraphQ[" H2 "], "
                   "ConnectedHypergraphQ[{{1,2},{2,3}}], ConnectedHypergraphQ[Hypergraph[{},{}]], "
                   "ConnectedHypergraphQ[x]}",
                   "{False, True, True, False, False}", 0);
    assert_eval_eq("HyperedgeConnectedComponents[" H1 "]", "{{1, 2, 3}, {4}}", 0);
    assert_eval_eq("HyperedgeConnectedComponents[" H2 ", 2]", "{{1, 2, 3}, {4}}", 0);
    /* Hyperedges smaller than s are in no s-walk. */
    assert_eval_eq("HyperedgeConnectedComponents[" H2 ", 3]", "{{1}, {2}, {3}}", 0);
    assert_eval_eq("HyperedgeConnectedComponents[Hypergraph[{{},{1}}]]", "{{2}}", 0);
    assert_eval_eq("{HyperedgeDistance[" H1 ", 1, 3], HyperedgeDistance[" H1 ", 1, 4], "
                   "HyperedgeDistance[" H1 ", 2, 2]}", "{2, Infinity, 0}", 0);
    assert_eval_eq("HyperedgeDistance[" H1 ", 1]", "{0, 1, 2, Infinity}", 0);
    assert_eval_eq("HyperedgeDistance[" H2 ", 1, All, 2]", "{0, 1, 2, Infinity}", 0);
    assert_eval_eq("HyperedgeDistance[" H2 ", 1, 4, 1]", "1", 0);
    assert_eval_eq("Head[HyperedgeDistance[" H1 ", 9, 1]]", "HyperedgeDistance", 0);
    assert_eval_eq("{HypergraphDistance[" H1 ", 1, 6], HypergraphDistance[" H1 ", 1, 7], "
                   "HypergraphDistance[" H1 ", 2, 2]}", "{3, Infinity, 0}", 0);
    assert_eval_eq("HypergraphDistance[" H1 ", 1]", "{0, 1, 1, 2, 3, 3, Infinity}", 0);
    /* Agrees with GraphDistance on the clique expansion. */
    assert_eval_eq("SeedRandom[5]; h = RandomHypergraph[{60, 40}, 3]; "
                   "HypergraphDistance[h, 1] === (GraphDistance[HypergraphCliqueExpansion[h], 1, #] & "
                   "/@ VertexList[h])", "True", 0);
    /* Line-graph BFS agrees with GraphDistance on the line graph, s = 1 and 2. */
    assert_eval_eq("SeedRandom[6]; h = RandomHypergraph[{30, 60}, 4]; "
                   "{HyperedgeDistance[h, 1] === (GraphDistance[HypergraphLineGraph[h], 1, #] & /@ Range[60]), "
                   "HyperedgeDistance[h, 1, All, 2] === (GraphDistance[HypergraphLineGraph[h, 2], 1, #] & /@ Range[60])}",
                   "{True, True}", 0);
    assert_eval_eq("SeedRandom[6]; h = RandomHypergraph[{30, 60}, 4]; "
                   "Length[HyperedgeConnectedComponents[h, 2]] === "
                   "Length[ConnectedComponents[HypergraphLineGraph[h, 2]]]", "True", 0);
}

static void test_edits(void) {
    assert_eval_eq("InputForm[HypergraphVertexAdd[" H1 ", {8, 1, 8}]]",
        "Hypergraph[{1, 2, 3, 4, 5, 6, 7, 8}, {{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]", 0);
    assert_eval_eq("VertexCount[HypergraphVertexAdd[" H1 ", x]]", "8", 0);
    assert_eval_eq("InputForm[HypergraphVertexDelete[" H1 ", 3]]",
        "Hypergraph[{1, 2, 4, 5, 6, 7}, {{4, 5, 6}, {7}}]", 0);
    assert_eval_eq("EdgeList[HypergraphVertexDelete[" H1 ", {4, 7}]]", "{{1, 2, 3}}", 0);
    assert_eval_eq("Head[HypergraphVertexDelete[" H1 ", 99]]", "HypergraphVertexDelete", 0);
    assert_eval_eq("InputForm[HypergraphEdgeAdd[" H1 ", {1, 10}]]",
        "Hypergraph[{1, 2, 3, 4, 5, 6, 7, 10}, {{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}, {1, 10}}]", 0);
    assert_eval_eq("EdgeList[HypergraphEdgeAdd[" H1 ", {{8}, {3, 4}}]]",
        "{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}, {8}, {3, 4}}", 0);
    assert_eval_eq("EdgeList[HypergraphEdgeDelete[Hypergraph[{{1,2},{2,3},{1,2}}], {1, 2}]]",
                   "{{2, 3}}", 0);
    assert_eval_eq("Head[HypergraphEdgeDelete[" H1 ", {4, 3}]]", "HypergraphEdgeDelete", 0);
    assert_eval_eq("InputForm[Subhypergraph[" H1 ", {7, 4, 3, 1, 2, 99}]]",
        "Hypergraph[{1, 2, 3, 4, 7}, {{1, 2, 3}, {3, 4}, {7}}]", 0);
    assert_eval_eq("InputForm[HypergraphRestriction[" H1 ", {1, 2, 3, 4}]]",
        "Hypergraph[{1, 2, 3, 4}, {{1, 2, 3}, {3, 4}, {4}}]", 0);
    assert_eval_eq("EdgeList[HypergraphRestriction[Hypergraph[{{1,2,1,3}}], {1, 3}]]",
                   "{{1, 1, 3}}", 0);
}

static void test_random(void) {
    assert_eval_eq("SeedRandom[1]; h = RandomHypergraph[{50, 200}, 4]; "
                   "{VertexCount[h], EdgeCount[h], UniformHypergraphQ[h, 4], "
                   "AllTrue[EdgeList[h], OrderedQ[#] && Length[Union[#]] == 4 &]}",
                   "{50, 200, True, True}", 0);
    assert_eval_eq("SeedRandom[9]; ra = RandomHypergraph[{20, 5}, 3]; SeedRandom[9]; "
                   "rb = RandomHypergraph[{20, 5}, 3]; ra === rb", "True", 0);
    assert_eval_eq("EdgeList[RandomHypergraph[{5, 3}, 5]]",
                   "{{1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}}", 0);
    assert_eval_eq("EdgeList[RandomHypergraph[{5, 2}, 0]]", "{{}, {}}", 0);
    assert_eval_eq("Length[RandomHypergraph[{10, 3}, 2, 4]]", "4", 0);
    assert_eval_eq("Head[RandomHypergraph[{3, 3}, 4]]", "RandomHypergraph", 0);
    /* FR form: arities and counts, vertices drawn with replacement. */
    assert_eval_eq("SeedRandom[2]; h = RandomHypergraph[{10, {{4, 3}, {2, 5}}}]; "
                   "{EdgeCount[h], HyperedgeSizes[h], VertexCount[h] <= 10}",
                   "{6, {3, 3, 3, 3, 5, 5}, True}", 0);
}

static void test_transversals(void) {
    assert_eval_eq("TransversalHypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]",
        "{{1, 4, 7}, {2, 4, 7}, {3, 4, 7}, {3, 5, 7}, {3, 6, 7}}", 0);
    assert_eval_eq("TransversalHypergraph[{{a,b},{b,c},{c,d},{d,a}}]", "{{a, c}, {b, d}}", 0);
    assert_eval_eq("InputForm[TransversalHypergraph[Hypergraph[{1,2,3},{{1,2},{2,3}}]]]",
                   "Hypergraph[{1, 2, 3}, {{2}, {1, 3}}]", 0);
    /* Degenerate inputs the FR function leaves unevaluated. */
    assert_eval_eq("{TransversalHypergraph[{}], TransversalHypergraph[{{}}], "
                   "TransversalHypergraph[{{1},{1,2}}], TransversalHypergraph[{{1,2},{1,2}}]}",
                   "{{{}}, {}, {{1}}, {{1}, {2}}}", 0);
    /* Tr(Tr(H)) = min(H) for a Sperner family. */
    assert_eval_eq("Sort /@ TransversalHypergraph[TransversalHypergraph[{{1,2},{2,3,4},{4,5}}]]",
                   "{{1, 2}, {4, 5}, {2, 3, 4}}", 0);
    /* Brute force on random instances (plain vs. repeated-vertex hyperedges). */
    assert_eval_eq(
        "iq[a_, b_] := Intersection[a, b] =!= {}; sq[a_, b_] := Complement[b, a] === {}; "
        "bf[e_, n_] := Module[{ts = Select[Subsets[Range[n]], Function[s, AllTrue[e, iq[#, s] &]]]}, "
        "  Select[ts, Function[s, NoneTrue[ts, Function[t, t =!= s && sq[s, t]]]]]]; "
        "SeedRandom[7]; "
        "{And @@ Table[With[{e = EdgeList[RandomHypergraph[{8, 6}, 3]]}, "
        "   Sort[Sort /@ TransversalHypergraph[e]] === Sort[bf[e, 8]]], {25}], "
        " And @@ Table[With[{e = EdgeList[RandomHypergraph[{9, {{5, 2}, {4, 3}}}]]}, "
        "   Sort[Sort /@ TransversalHypergraph[e]] === Sort[bf[e, 9]]], {25}], "
        " And @@ Table[With[{e = EdgeList[RandomHypergraph[{10, 8}, 3]]}, "
        "   With[{t = FindMinimumTransversal[e]}, AllTrue[e, iq[#, t] &] && "
        "     Length[t] === Min[Length /@ bf[e, 10]]]], {25}]}",
        "{True, True, True}", 0);
    assert_eval_eq("FindMinimumTransversal[{{1,2,3},{3,4},{4,5,6},{7}}]", "{3, 4, 7}", 0);
    assert_eval_eq("FindMinimumTransversal[{}]", "{}", 0);
    assert_eval_eq("Head[FindMinimumTransversal[{{1},{}}]]", "FindMinimumTransversal", 0);
    /* Minimum transversal of the 5-cycle's edges is 3; of disjoint pairs, one each. */
    assert_eval_eq("Length[FindMinimumTransversal[EdgeList[Hypergraph[CycleGraph[5]]]]]", "3", 0);
    assert_eval_eq("Length[FindMinimumTransversal[Table[{2 i - 1, 2 i}, {i, 500}]]]", "500", 0);
}

static void test_memo(void) {
    /* In-place Part assignment must not be answered from a stale memo entry. */
    assert_eval_eq("hm = " H1 "; {VertexDegree[hm, 3], EdgeCount[hm]}", "{2, 4}", 0);
    assert_eval_eq("hm[[2]] = {{1, 2}}; {VertexDegree[hm, 3], EdgeCount[hm], HypergraphQ[hm]}",
                   "{0, 1, True}", 0);
    assert_eval_eq("hm[[2]] = {{1, 99}}; HypergraphQ[hm]", "False", 0);
    /* More live hypergraphs than memo slots. */
    assert_eval_eq("hs = Table[Hypergraph[{Range[k]}], {k, 1, 9}]; VertexCount /@ hs",
                   "{1, 2, 3, 4, 5, 6, 7, 8, 9}", 0);
    assert_eval_eq("VertexDegree /@ Take[hs, 3]", "{{1}, {1, 1}, {1, 1, 1}}", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_construction);
    TEST(test_accessors);
    TEST(test_sizes);
    TEST(test_dual_and_expansions);
    TEST(test_line_graph);
    TEST(test_components_and_distances);
    TEST(test_edits);
    TEST(test_random);
    TEST(test_transversals);
    TEST(test_memo);

    printf("All hypergraph tests passed!\n");
    return 0;
}
