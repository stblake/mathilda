/* test_graph_metrics.c - distance, centrality, clustering and generator
 * builtins (src/graph/gmet_*.c).
 *
 * Every expected value below was checked against Mathematica 15 with
 * wolframscript on the same input (machine reals compared as Round[10^6 x]);
 * the few intentional deviations are marked: KirchhoffMatrix is dense (Wolfram
 * returns a SparseArray), BetweennessCentrality leaves a mixed graph
 * unevaluated, and the WheelGraph[3] / PetersenGraph[5, 5] multigraphs are left
 * unevaluated because Graph is simple.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "print.h"
#include "graph.h"
#include "graph_metrics.h"
#include "test_utils.h"
#include <stdlib.h>

static void test_distances(void) {
    assert_eval_eq("GraphDistanceMatrix[PathGraph[Range[4]]]",
                   "{{0, 1, 2, 3}, {1, 0, 1, 2}, {2, 1, 0, 1}, {3, 2, 1, 0}}", 0);
    assert_eval_eq("GraphDistanceMatrix[Graph[{1->2, 2->3, 3->1, 3->4}]]",
                   "{{0, 1, 2, 3}, {2, 0, 1, 2}, {1, 2, 0, 1}, {Infinity, Infinity, Infinity, 0}}", 0);
    assert_eval_eq("GraphDistanceMatrix[Graph[{1, 2, 3}, {1<->2, 2<->3, 1<->3}, EdgeWeight->{1, 5, 2}]]",
                   "{{0.0, 1.0, 2.0}, {1.0, 0.0, 3.0}, {2.0, 3.0, 0.0}}", 0);
    assert_eval_eq("GraphDistanceMatrix[Graph[{1, 2, 3}, {1<->2, 2<->3}, EdgeWeight->{1.5, 2}], 2]",
                   "{{0.0, 1.5, Infinity}, {1.5, 0.0, 2.0}, {Infinity, 2.0, 0.0}}", 0);
    assert_eval_eq("GraphDistanceMatrix[PathGraph[Range[5]], 2]",
                   "{{0, 1, 2, Infinity, Infinity}, {1, 0, 1, 2, Infinity}, {2, 1, 0, 1, 2}, {Infinity, 2, 1, 0, 1}, {Infinity, Infinity, 2, 1, 0}}", 0);
    assert_eval_eq("NDArrayQ[GraphDistanceMatrix[PathGraph[Range[4]]]]",
                   "True", 0);
    assert_eval_eq("NDArrayQ[GraphDistanceMatrix[Graph[{1, 2, 3}, {1<->2, 2<->3, 1<->3}, EdgeWeight->{1, 5, 2}]]]",
                   "True", 0);
    assert_eval_eq("GraphDistanceMatrix[Graph[{1, 2, 3}, {1<->2, 2<->3}, EdgeWeight->{x, 2}]]",
                   "GraphDistanceMatrix[Graph[<3 vertices, 2 edges>]]", 0);
    assert_eval_eq("GraphDistance[Graph[{1<->2, 2<->3, 3<->4, 4<->1, 1<->3, 5<->6}], 1]",
                   "{0, 1, 1, 1, Infinity, Infinity}", 0);
    assert_eval_eq("GraphDistance[Graph[{1->2, 2->3, 3->1, 3->4}], 1]",
                   "{0, 1, 2, 3}", 0);
    assert_eval_eq("GraphDistance[Graph[{1, 2, 3}, {1<->2, 2<->3, 1<->3}, EdgeWeight->{1, 5, 2}], 2]",
                   "{1.0, 0, 3.0}", 0);
    assert_eval_eq("GraphDistance[PathGraph[Range[5]], 1, 4]",
                   "3", 0);
    assert_eval_eq("NDArrayQ[GraphDistance[PathGraph[Range[10]], 1]]",
                   "True", 0);
    assert_eval_eq("{VertexEccentricity[Graph[{1<->2, 2<->3, 3<->4, 4<->1, 1<->3, 5<->6}], 1], VertexEccentricity[PathGraph[Range[5]], 1], VertexEccentricity[Graph[{1->2, 2->3, 3->1, 3->4}], 4]}",
                   "{1, 4, 0}", 0);
    assert_eval_eq("VertexEccentricity[Graph[{1, 2, 3}, {1<->2}, EdgeWeight->{2}], 1]",
                   "Infinity", 0);
    assert_eval_eq("{GraphDiameter[PathGraph[Range[5]]], GraphRadius[PathGraph[Range[5]]], GraphCenter[PathGraph[Range[5]]], GraphPeriphery[PathGraph[Range[5]]], MeanGraphDistance[PathGraph[Range[5]]]}",
                   "{4, 2, {3}, {1, 5}, 2}", 0);
    assert_eval_eq("{GraphDiameter[CycleGraph[5]], MeanGraphDistance[Graph[{1->2, 2->3, 3->1}]], MeanGraphDistance[Graph[{1<->2, 3<->4}]]}",
                   "{2, 3/2, Infinity}", 0);
    assert_eval_eq("{GraphDiameter[Graph[{1->2, 2->3, 3->1, 3->4}]], GraphRadius[Graph[{1->2, 2->3, 3->1, 3->4}]], GraphCenter[Graph[{1->2, 2->3, 3->1, 3->4}]], GraphPeriphery[Graph[{1->2, 2->3, 3->1, 3->4}]]}",
                   "{Infinity, Infinity, {}, {}}", 0);
    assert_eval_eq("{GraphDiameter[Graph[{1, 2, 3, 4, 5}, {1->2, 2->3, 3->1, 3->4, 4->5, 1->5}, EdgeWeight->{2, 1, 3, 1, 1, 7}]], GraphRadius[Graph[{1, 2, 3, 4, 5}, {1->2, 2->3, 3->1, 3->4, 4->5, 1->5}, EdgeWeight->{2, 1, 3, 1, 1, 7}]], GraphCenter[Graph[{1, 2, 3, 4, 5}, {1->2, 2->3, 3->1, 3->4, 4->5, 1->5}, EdgeWeight->{2, 1, 3, 1, 1, 7}]], GraphPeriphery[Graph[{1, 2, 3, 4, 5}, {1->2, 2->3, 3->1, 3->4, 4->5, 1->5}, EdgeWeight->{2, 1, 3, 1, 1, 7}]]}",
                   "{Infinity, 4.0, {2}, {4, 5}}", 0);
    assert_eval_eq("{GraphDiameter[Graph[{1}, {}]], GraphCenter[Graph[{1}, {}]], MeanGraphDistance[Graph[{1}, {}]]}",
                   "{0, {1}, MeanGraphDistance[Graph[<1 vertex, 0 edges>]]}", 0);
    assert_eval_eq("{GraphDensity[PathGraph[Range[5]]], GraphDensity[Graph[{1->2, 2->3, 3->1, 3->4}]], GraphDensity[Graph[{1, 2, 3}, {1->2, 2<->3}]], GraphDensity[Graph[{1}, {}]]}",
                   "{2/5, 1/3, 1/2, GraphDensity[Graph[<1 vertex, 0 edges>]]}", 0);
}

static void test_density_kirchhoff(void) {
    assert_eval_eq("KirchhoffMatrix[Graph[{1->2, 2->3, 3->1, 3->4}]]",
                   "{{2, -1, 0, 0}, {0, 2, -1, 0}, {-1, 0, 3, -1}, {0, 0, 0, 1}}", 0);
    assert_eval_eq("KirchhoffMatrix[Graph[{1->2, 2<->3}]]",
                   "{{1, -1, 0}, {0, 2, -1}, {0, -1, 1}}", 0);
    assert_eval_eq("NDArrayQ[KirchhoffMatrix[CycleGraph[5]]]",
                   "True", 0);
}

static void test_path_centralities(void) {
    assert_eval_eq("{DegreeCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]], DegreeCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], \"In\"], DegreeCentrality[Graph[{1->2, 2->3, 3->1, 3->4}], \"Out\"], DegreeCentrality[Graph[{1->2, 2<->3}]]}",
                   "{{2, 2, 3, 1}, {1, 1, 1, 1}, {1, 1, 2, 0}, {1, 3, 2}}", 0);
    assert_eval_eq("NDArrayQ[DegreeCentrality[CycleGraph[6]]]",
                   "True", 0);
    assert_eval_eq("Round[10^6 ClosenessCentrality[Graph[{1<->2, 2<->3, 3<->4, 4<->1, 1<->3, 5<->6}]]]",
                   "{1000000, 750000, 1000000, 750000, 1000000, 1000000}", 0);
    assert_eval_eq("Round[10^6 ClosenessCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]]",
                   "{500000, 600000, 750000, 0}", 0);
    assert_eval_eq("Round[10^6 ClosenessCentrality[Graph[{1, 2, 3}, {1<->2, 2<->3, 1<->3}, EdgeWeight->{1.5, 5, 2}]]]",
                   "{571429, 400000, 363636}", 0);
    assert_eval_eq("Round[10^6 EccentricityCentrality[Graph[{1, 2, 3, 4, 5}, {1->2, 2->3, 3->1, 3->4, 4->5, 1->5}, EdgeWeight->{2, 1, 3, 1, 1, 7}]]]",
                   "{200000, 250000, 200000, 1000000, 0}", 0);
    assert_eval_eq("NDArrayQ[ClosenessCentrality[PathGraph[Range[10]]]]",
                   "True", 0);
    assert_eval_eq("BetweennessCentrality[PathGraph[Range[5]]]",
                   "{0.0, 3.0, 4.0, 3.0, 0.0}", 0);
    assert_eval_eq("BetweennessCentrality[Graph[{1->2, 2->3, 3->1, 3->4}]]",
                   "{1.0, 2.0, 3.0, 0.0}", 0);
    assert_eval_eq("BetweennessCentrality[Graph[{1, 2, 3, 4}, {1<->2, 2<->3, 3<->4, 4<->1}, EdgeWeight->{1, 1, 1, 5}]]",
                   "{0.5, 0.5, 0.5, 0.5}", 0);
    assert_eval_eq("BetweennessCentrality[Graph[{1<->2, 2->3}]]",
                   "BetweennessCentrality[Graph[<3 vertices, 2 edges>]]", 0);
    assert_eval_eq("{BetweennessCentrality[StarGraph[5]], BetweennessCentrality[Graph[{1}, {}]]}",
                   "{{6.0, 0.0, 0.0, 0.0, 0.0}, {0.0}}", 0);
    assert_eval_eq("EdgeBetweennessCentrality[Graph[{1<->2, 2<->3, 3<->4, 4<->1, 1<->3, 5<->6}]]",
                   "{3.0, 3.0, 3.0, 3.0, 2.0, 2.0}", 0);
    assert_eval_eq("EdgeBetweennessCentrality[Graph[{1, 2, 3, 4}, {1<->2, 2<->3, 3<->4, 4<->1}, EdgeWeight->{1, 1, 1, 3}]]",
                   "{5.0, 7.0, 5.0, 1.0}", 0);
    assert_eval_eq("EdgeBetweennessCentrality[Graph[{1<->2, 2<->3, 3->4}]]",
                   "{5.0, 6.0, 3.0}", 0);
    assert_eval_eq("Round[10^6 PageRankCentrality[Graph[{1->2, 1->3, 2->3}], 0.85]]",
                   "{197580, 281551, 520869}", 0);
}

static void test_spectral_centralities(void) {
    assert_eval_eq("Round[10^6 PageRankCentrality[Graph[{1, 2, 3}, {1->2}]]]",
                   "{259740, 480519, 259740}", 0);
    assert_eval_eq("Round[10^6 PageRankCentrality[PathGraph[Range[5]], 1/2]]",
                   "{158333, 233333, 216667, 233333, 158333}", 0);
    assert_eval_eq("Round[10^6 EigenvectorCentrality[Graph[{1<->2, 2<->3, 3<->4, 4<->1, 1<->3, 2<->4, 5<->6, 6<->7, 7<->5, 8<->9}]]]",
                   "{125000, 125000, 125000, 125000, 111111, 111111, 111111, 83333, 83333}", 0);
    assert_eval_eq("Round[10^6 EigenvectorCentrality[Graph[{1->2, 2->3, 3->4, 4->1, 1->3}]]]",
                   "{220744, 180827, 328956, 269472}", 0);
    assert_eval_eq("Round[10^6 EigenvectorCentrality[Graph[{1->2, 2->3, 3->4, 4->1, 1->3}], \"Out\"]]",
                   "{328956, 180827, 220744, 269472}", 0);
    assert_eval_eq("Round[10^6 EigenvectorCentrality[Graph[{1->2, 2->1, 2->3, 3->4, 4->5, 5->3}]]]",
                   "{166667, 166667, 222222, 222222, 222222}", 0);
    assert_eval_eq("EigenvectorCentrality[Graph[{1->2, 1->3, 2->3}]]",
                   "{0.0, 0.0, 0.0}", 0);
    assert_eval_eq("Round[10^6 KatzCentrality[PathGraph[Range[5]], 0.1]]",
                   "{1123711, 1237113, 1247423, 1237113, 1123711}", 0);
    assert_eval_eq("Round[10^6 KatzCentrality[Graph[{1->2, 2->3}], 1/2]]",
                   "{1000000, 1500000, 1750000}", 0);
    assert_eval_eq("Round[10^6 KatzCentrality[CompleteGraph[4], 0.5]]",
                   "{-2000000, -2000000, -2000000, -2000000}", 0);
    /* singular up to rounding: unevaluated (Mathematica returns ~1.8*10^16) */
    assert_eval_eq("Head[KatzCentrality[CompleteGraph[4], 1/3]]", "KatzCentrality", 0);
    assert_eval_eq("Head[KatzCentrality[CycleGraph[4], 0.5]]", "KatzCentrality", 0);
    assert_eval_eq("Round[KatzCentrality[CompleteGraph[4], 0.3333]]", "{10000, 10000, 10000, 10000}", 0);
    assert_eval_eq("Round[10^6 KatzCentrality[PathGraph[Range[4]], 0.2, {1, 2, 3, 4}]]",
                   "{1651543, 3257713, 4637024, 4927405}", 0);
    assert_eval_eq("KatzCentrality[CycleGraph[4], 0]",
                   "{1, 1, 1, 1}", 0);
    assert_eval_eq("KatzCentrality[CompleteGraph[3], 0.5]",
                   "KatzCentrality[Graph[<3 vertices, 3 edges>], 0.5]", 0);
    assert_eval_eq("Round[10^6 HITSCentrality[Graph[{1->2, 1->3, 2->3}]]]",
                   "{{0, 381966, 618034}, {1000000, 618034, 0}}", 0);
    assert_eval_eq("Round[10^6 HITSCentrality[PathGraph[Range[5]]]]",
                   "{{166667, 166667, 333333, 166667, 166667}, {166667, 500000, 333333, 500000, 166667}}", 0);
    assert_eval_eq("Round[10^6 HITSCentrality[Graph[{1, 2, 3, 4, 5}, {1->2, 2->3, 3->1, 3->4, 4->5, 1->5}]]]",
                   "{{250000, 190983, 0, 250000, 309017}, {500000, 0, 500000, 309017, 0}}", 0);
    assert_eval_eq("HITSCentrality[Graph[{1->2, 2<->3}]]",
                   "{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}", 0);
}

static void test_clustering(void) {
    assert_eval_eq("{GraphTriangleCount[CompleteGraph[30]], GraphTriangleCount[Graph[{1->2, 2->3, 3->1, 3->4}]], GraphTriangleCount[Graph[{1->2, 2->1, 2->3, 3->2, 3->1, 1->3}]], GraphTriangleCount[Graph[{1->2, 2->3, 1->3}]]}",
                   "{4060, 1, 2, 0}", 0);
    assert_eval_eq("GraphTriangleCount[Graph[{1<->2, 2->3, 3->1}]]",
                   "GraphTriangleCount[Graph[<3 vertices, 3 edges>]]", 0);
    assert_eval_eq("{LocalClusteringCoefficient[Graph[{1<->2, 2<->3, 3<->4, 4<->1, 1<->3, 5<->6}]], LocalClusteringCoefficient[Graph[{1->2, 2->3, 3->1, 3->4}]], LocalClusteringCoefficient[Graph[{1<->2, 2<->3, 3<->1, 3<->4}], 3]}",
                   "{{2/3, 1, 2/3, 1, 0, 0}, {1, 1, 1/2, 0}, 1/3}", 0);
    assert_eval_eq("{GlobalClusteringCoefficient[Graph[{1<->2, 2<->3, 3<->4, 4<->1, 1<->3, 5<->6}]], GlobalClusteringCoefficient[Graph[{1->2, 2->3, 3->1, 3->4}]], GlobalClusteringCoefficient[PathGraph[Range[4]]]}",
                   "{3/4, 3/4, 0}", 0);
    assert_eval_eq("{MeanClusteringCoefficient[Graph[{1<->2, 2<->3, 3<->4, 4<->1, 1<->3, 5<->6}]], MeanClusteringCoefficient[Graph[{1<->2, 2<->3, 3<->1, 3<->4, 4<->5, 5<->3, 5<->6}]], MeanClusteringCoefficient[Graph[{}, {}]]}",
                   "{5/9, 11/18, 0}", 0);
}

static void test_generators(void) {
    assert_eval_eq("EdgeList[WheelGraph[5]]",
                   "{1 <-> 2, 1 <-> 3, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 5, 3 <-> 4, 4 <-> 5}", 0);
    assert_eval_eq("EdgeList[HypercubeGraph[3]]",
                   "{1 <-> 2, 1 <-> 3, 1 <-> 5, 2 <-> 4, 2 <-> 6, 3 <-> 4, 3 <-> 7, 4 <-> 8, 5 <-> 6, 5 <-> 7, 6 <-> 8, 7 <-> 8}", 0);
    assert_eval_eq("EdgeList[GridGraph[{2, 3}]]",
                   "{1 <-> 2, 1 <-> 3, 2 <-> 4, 3 <-> 4, 3 <-> 5, 4 <-> 6, 5 <-> 6}", 0);
    assert_eval_eq("EdgeList[GridGraph[{2, 2, 2}]]",
                   "{1 <-> 2, 1 <-> 3, 1 <-> 5, 2 <-> 4, 2 <-> 6, 3 <-> 4, 3 <-> 7, 4 <-> 8, 5 <-> 6, 5 <-> 7, 6 <-> 8, 7 <-> 8}", 0);
    assert_eval_eq("{EdgeList[KaryTree[7, 3]], EdgeList[CompleteKaryTree[3, 3]]}",
                   "{{1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 5, 2 <-> 6, 2 <-> 7}, {1 <-> 2, 1 <-> 3, 1 <-> 4, 2 <-> 5, 2 <-> 6, 2 <-> 7, 3 <-> 8, 3 <-> 9, 3 <-> 10, 4 <-> 11, 4 <-> 12, 4 <-> 13}}", 0);
    assert_eval_eq("{EdgeList[CirculantGraph[6, {1, 2}]], EdgeList[CirculantGraph[5, {7}]]}",
                   "{{1 <-> 2, 1 <-> 3, 1 <-> 5, 1 <-> 6, 2 <-> 3, 2 <-> 4, 2 <-> 6, 3 <-> 4, 3 <-> 5, 4 <-> 5, 4 <-> 6, 5 <-> 6}, {1 <-> 3, 1 <-> 4, 2 <-> 4, 2 <-> 5, 3 <-> 5}}", 0);
    assert_eval_eq("EdgeList[PetersenGraph[]]",
                   "{1 <-> 3, 1 <-> 4, 1 <-> 6, 2 <-> 4, 2 <-> 5, 2 <-> 7, 3 <-> 5, 3 <-> 8, 4 <-> 9, 5 <-> 10, 6 <-> 7, 6 <-> 10, 7 <-> 8, 8 <-> 9, 9 <-> 10}", 0);
    assert_eval_eq("EdgeList[PetersenGraph[4, 2]]",
                   "{1 <-> 3, 1 <-> 5, 2 <-> 4, 2 <-> 6, 3 <-> 7, 4 <-> 8, 5 <-> 6, 5 <-> 8, 6 <-> 7, 7 <-> 8}", 0);
    assert_eval_eq("EdgeList[TuranGraph[7, 3]]",
                   "{1 <-> 4, 1 <-> 5, 1 <-> 6, 1 <-> 7, 2 <-> 4, 2 <-> 5, 2 <-> 6, 2 <-> 7, 3 <-> 4, 3 <-> 5, 3 <-> 6, 3 <-> 7, 4 <-> 6, 4 <-> 7, 5 <-> 6, 5 <-> 7}", 0);
    assert_eval_eq("{EdgeList[CompleteGraph[{2, 3}]], EdgeList[CompleteGraph[{3}]], VertexCount[CompleteGraph[4]]}",
                   "{{1 <-> 3, 1 <-> 4, 1 <-> 5, 2 <-> 3, 2 <-> 4, 2 <-> 5}, {1 <-> 2, 1 <-> 3, 2 <-> 3}, 4}", 0);
    assert_eval_eq("{EdgeList[HararyGraph[3, 7]], EdgeList[HararyGraph[5, 8]]}",
                   "{{1 <-> 2, 1 <-> 4, 1 <-> 5, 1 <-> 7, 2 <-> 3, 2 <-> 6, 3 <-> 4, 3 <-> 7, 4 <-> 5, 5 <-> 6, 6 <-> 7}, {1 <-> 2, 1 <-> 3, 1 <-> 5, 1 <-> 7, 1 <-> 8, 2 <-> 3, 2 <-> 4, 2 <-> 6, 2 <-> 8, 3 <-> 4, 3 <-> 5, 3 <-> 7, 4 <-> 5, 4 <-> 6, 4 <-> 8, 5 <-> 6, 5 <-> 7, 6 <-> 7, 6 <-> 8, 7 <-> 8}}", 0);
    assert_eval_eq("{WheelGraph[3], PetersenGraph[5, 5], HararyGraph[1, 5], GridGraph[{0, 3}], CompleteGraph[{2, 0}]}",
                   "{WheelGraph[3], PetersenGraph[5, 5], HararyGraph[1, 5], GridGraph[{0, 3}], CompleteGraph[{2, 0}]}", 0);
}

/* The per-graph caches (result cache, distance summary) hold a reference to
 * the graph node, so an in-place Part assignment must unshare it and the next
 * call must see the new graph rather than a stale cached answer. */
static void test_cache_safety(void) {
    assert_eval_eq("gm = PathGraph[Range[4]]; BetweennessCentrality[gm]",
                   "{0.0, 2.0, 2.0, 0.0}", 0);
    assert_eval_eq("gm[[2]] = {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}; BetweennessCentrality[gm]",
                   "{0.5, 0.5, 0.5, 0.5}", 0);
    assert_eval_eq("gm = PathGraph[Range[4]]; {GraphDiameter[gm], ClosenessCentrality[gm][[1]]}",
                   "{3, 0.5}", 0);
    assert_eval_eq("gm[[2]] = {1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}; "
                   "{GraphDiameter[gm], ClosenessCentrality[gm][[1]]}",
                   "{2, 0.75}", 0);
    /* More live graphs than cache slots: every answer still correct. */
    assert_eval_eq("Table[GraphDiameter[PathGraph[Range[k]]], {k, 2, 40}] == Range[1, 39]",
                   "True", 0);
    assert_eval_eq("Table[Total[DegreeCentrality[CycleGraph[k]]], {k, 3, 30}] == 2 Range[3, 30]",
                   "True", 0);
}

/* Larger graphs: the threaded and bit-parallel paths agree with simple
 * identities (sum of distances, handshake-style totals). */
static void test_scale_identities(void) {
    /* Grid distances: Manhattan metric. */
    assert_eval_eq("g = GridGraph[{30, 20}]; d = GraphDistanceMatrix[g]; "
                   "{Dimensions[d], NDArrayQ[d], d[[1, 600]], Max[d], GraphDiameter[g]}",
                   "{{600, 600}, True, 48, 48, 48}", 0);
    /* Betweenness sums to sum over pairs of (distance - 1). */
    assert_eval_eq("g = GridGraph[{12, 9}]; d = GraphDistanceMatrix[g]; "
                   "Round[Total[BetweennessCentrality[g]] - (Total[Flatten[d]] - 108*107)/2]",
                   "0", 0);
    /* Closeness from the matrix. */
    assert_eval_eq("g = CirculantGraph[500, {1, 7, 31}]; d = GraphDistanceMatrix[g]; "
                   "Max[Abs[ClosenessCentrality[g] - 499/Total[d, {2}]]] < 10^-12",
                   "True", 0);
    /* Triangle counts: K_n has Binomial[n, 3]; directed K_n has 2 Binomial[n, 3]. */
    assert_eval_eq("{GraphTriangleCount[CompleteGraph[60]], "
                   "GraphTriangleCount[Graph[Range[12], Flatten[Table[If[i != j, i -> j, Nothing], {i, 12}, {j, 12}]]]]}",
                   "{34220, 440}", 0);
    /* PageRank sums to 1; eigenvector centrality of a regular graph is uniform. */
    assert_eval_eq("{Round[10^9 Total[PageRankCentrality[HypercubeGraph[10]]]], "
                   "Round[10^9 Max[EigenvectorCentrality[HypercubeGraph[10]]] 1024]}",
                   "{1000000000, 1000000000}", 0);
    /* Weighted all-pairs on a path: distances are prefix sums. */
    assert_eval_eq("g = Graph[Range[300], Table[i <-> i + 1, {i, 299}], EdgeWeight -> Table[1/2, {299}]]; "
                   "{GraphDistanceMatrix[g][[1, 300]], GraphDiameter[g], MeanGraphDistance[g]}",
                   "{149.5, 149.5, 50.1667}", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_distances);
    TEST(test_density_kirchhoff);
    TEST(test_path_centralities);
    TEST(test_spectral_centralities);
    TEST(test_clustering);
    TEST(test_generators);
    TEST(test_cache_safety);
    TEST(test_scale_identities);

    printf("All graph metrics tests passed!\n");
    return 0;
}
