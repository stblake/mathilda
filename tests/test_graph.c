/* test_graph.c - Phase 1 graph subsystem tests.
 *
 * Covers construction/normalization (all four edge sugars), vertex derivation,
 * directed-by-default, rejection of self-loops / parallel edges / 3-arg edges /
 * unknown endpoints, the GraphQ predicate, terse-summary printing, and the
 * InputForm round-trip through the parser.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "print.h"
#include "test_utils.h"
#include <stdlib.h>

/* ---- Normalization + FullForm (edge sugar -> canonical edges) ------------- */
static void test_edge_sugar_normalization(void) {
    /* All four accepted edge forms canonicalize to Directed/UndirectedEdge. */
    assert_eval_eq("Graph[{1,2},{1->2}]",
                   "Graph[List[1, 2], List[DirectedEdge[1, 2]]]", 1);
    assert_eval_eq("Graph[{1,2},{DirectedEdge[1,2]}]",
                   "Graph[List[1, 2], List[DirectedEdge[1, 2]]]", 1);
    assert_eval_eq("Graph[{1,2},{1<->2}]",
                   "Graph[List[1, 2], List[UndirectedEdge[1, 2]]]", 1);
    assert_eval_eq("Graph[{1,2},{TwoWayRule[1,2]}]",
                   "Graph[List[1, 2], List[UndirectedEdge[1, 2]]]", 1);
    assert_eval_eq("Graph[{1,2},{UndirectedEdge[1,2]}]",
                   "Graph[List[1, 2], List[UndirectedEdge[1, 2]]]", 1);
}

/* ---- Vertex derivation (Graph[edges], directed by default) ---------------- */
static void test_vertex_derivation(void) {
    /* Vertices derived in first-appearance order; Rule defaults to directed. */
    assert_eval_eq("Graph[{1->2,2->3,3->1}]",
                   "Graph[List[1, 2, 3], "
                   "List[DirectedEdge[1, 2], DirectedEdge[2, 3], DirectedEdge[3, 1]]]", 1);
    /* Derivation preserves the order endpoints first appear, not sorted. */
    assert_eval_eq("Graph[{3->1,1->2}]",
                   "Graph[List[3, 1, 2], "
                   "List[DirectedEdge[3, 1], DirectedEdge[1, 2]]]", 1);
}

/* ---- Terse summary printing (standard form) ------------------------------- */
static void test_summary_printing(void) {
    assert_eval_eq("Graph[{1,2,3},{1->2,2->3}]", "Graph[<3 vertices, 2 edges>]", 0);
    assert_eval_eq("Graph[{1,2},{1->2}]",        "Graph[<2 vertices, 1 edge>]", 0);
    assert_eval_eq("Graph[{1},{}]",              "Graph[<1 vertex, 0 edges>]", 0);
}

/* ---- GraphQ truth table --------------------------------------------------- */
static void test_graphq(void) {
    assert_eval_eq("GraphQ[Graph[{1,2},{1->2}]]", "True", 0);
    assert_eval_eq("GraphQ[Graph[{1,2},{1<->2}]]", "True", 0);
    assert_eval_eq("GraphQ[Graph[{1,2,3},{1->2,2->3}]]", "True", 0);
    /* Not a graph at all. */
    assert_eval_eq("GraphQ[5]", "False", 0);
    assert_eval_eq("GraphQ[foo]", "False", 0);
}

/* ---- Rejection of malformed graphs (stay unevaluated -> GraphQ False) ------ */
static void test_rejections(void) {
    /* Self-loop. */
    assert_eval_eq("GraphQ[Graph[{1},{1->1}]]", "False", 0);
    /* Parallel/duplicate edges. */
    assert_eval_eq("GraphQ[Graph[{1,2},{1->2,1->2}]]", "False", 0);
    assert_eval_eq("GraphQ[Graph[{1,2},{1<->2,2<->1}]]", "False", 0);
    /* 3-argument edge (reserved for future edge tags). */
    assert_eval_eq("GraphQ[Graph[{1,2},{DirectedEdge[1,2,x]}]]", "False", 0);
    /* Edge endpoint absent from an explicit vertex list. */
    assert_eval_eq("GraphQ[Graph[{1,2},{1->3}]]", "False", 0);
    /* Anti-parallel *directed* edges are allowed (distinct). */
    assert_eval_eq("GraphQ[Graph[{1,2},{1->2,2->1}]]", "True", 0);
}

/* ---- InputForm round-trips through the parser ----------------------------- */
static void test_inputform_roundtrip(void) {
    static const char* inputs[] = {
        "Graph[{1,2},{1<->2}]",
        "Graph[{1,2,3},{1->2,2->3,3->1}]",
        "Graph[{a,b,c},{a<->b,b->c}]",
    };
    for (size_t i = 0; i < sizeof(inputs) / sizeof(inputs[0]); i++) {
        Expr* g = evaluate(parse_expression(inputs[i]));
        ASSERT(g != NULL);
        /* Wrap in InputForm and print: yields the literal constructor. */
        Expr* wrap_args[1] = { expr_copy(g) };
        Expr* wrap = expr_new_function(expr_new_symbol("InputForm"), wrap_args, 1);
        char* s = expr_to_string(wrap);
        /* Re-parse and evaluate: must reproduce an equal graph. */
        Expr* g2 = evaluate(parse_expression(s));
        ASSERT(expr_eq(g, g2));
        free(s);
        expr_free(wrap);
        expr_free(g);
        expr_free(g2);
    }
    /* Spot-check the exact InputForm text for the undirected case. */
    assert_eval_eq("InputForm[Graph[{1,2},{1<->2}]]", "Graph[{1, 2}, {1 <-> 2}]", 0);
    assert_eval_eq("InputForm[Graph[{1,2},{1->2}]]",  "Graph[{1, 2}, {1 -> 2}]", 0);
}

/* ---- Phase 2: query / representation builtins ----------------------------- */
static void test_query_builtins(void) {
    /* g = 1->2->3->4->1 (directed cycle). */
    const char* g = "Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]";
    char buf[256];

    snprintf(buf, sizeof(buf), "VertexList[%s]", g);
    assert_eval_eq(buf, "{1, 2, 3, 4}", 0);
    snprintf(buf, sizeof(buf), "EdgeList[%s]", g);
    assert_eval_eq(buf, "{1 -> 2, 2 -> 3, 3 -> 4, 4 -> 1}", 0);
    snprintf(buf, sizeof(buf), "VertexCount[%s]", g);
    assert_eval_eq(buf, "4", 0);
    snprintf(buf, sizeof(buf), "EdgeCount[%s]", g);
    assert_eval_eq(buf, "4", 0);

    /* Directed cycle: every vertex has in=out=1, total degree 2. */
    snprintf(buf, sizeof(buf), "VertexDegree[%s, 1]", g);
    assert_eval_eq(buf, "2", 0);
    snprintf(buf, sizeof(buf), "VertexInDegree[%s, 1]", g);
    assert_eval_eq(buf, "1", 0);
    snprintf(buf, sizeof(buf), "VertexOutDegree[%s, 1]", g);
    assert_eval_eq(buf, "1", 0);
    snprintf(buf, sizeof(buf), "VertexDegree[%s]", g);
    assert_eval_eq(buf, "{2, 2, 2, 2}", 0);

    /* Successor adjacency for a directed graph. */
    snprintf(buf, sizeof(buf), "AdjacencyList[%s, 1]", g);
    assert_eval_eq(buf, "{2}", 0);
    snprintf(buf, sizeof(buf), "AdjacencyList[%s]", g);
    assert_eval_eq(buf, "{{2}, {3}, {4}, {1}}", 0);

    snprintf(buf, sizeof(buf), "DirectedGraphQ[%s]", g);
    assert_eval_eq(buf, "True", 0);
}

static void test_query_undirected(void) {
    /* Path 1 <-> 2 <-> 3 (undirected). */
    const char* g = "Graph[{1,2,3},{1<->2,2<->3}]";
    char buf[256];

    /* Undirected: middle vertex has degree 2, ends degree 1. */
    snprintf(buf, sizeof(buf), "VertexDegree[%s]", g);
    assert_eval_eq(buf, "{1, 2, 1}", 0);
    /* Undirected neighbors go both ways. */
    snprintf(buf, sizeof(buf), "AdjacencyList[%s, 2]", g);
    assert_eval_eq(buf, "{1, 3}", 0);
    snprintf(buf, sizeof(buf), "DirectedGraphQ[%s]", g);
    assert_eval_eq(buf, "False", 0);
    /* In/out degree equal the degree for undirected graphs. */
    snprintf(buf, sizeof(buf), "VertexInDegree[%s, 2]", g);
    assert_eval_eq(buf, "2", 0);
    snprintf(buf, sizeof(buf), "VertexOutDegree[%s, 2]", g);
    assert_eval_eq(buf, "2", 0);
}

/* ---- Phase 3: matrix views ------------------------------------------------ */
static void test_matrix_views(void) {
    /* Directed 4-cycle: circulant adjacency matrix. */
    const char* dg = "Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]";
    char buf[256];
    snprintf(buf, sizeof(buf), "AdjacencyMatrix[%s]", dg);
    assert_eval_eq(buf,
        "{{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}", 0);

    /* Undirected edge -> symmetric matrix. */
    assert_eval_eq("AdjacencyMatrix[Graph[{1,2},{1<->2}]]",
                   "{{0, 1}, {1, 0}}", 0);

    /* Feeds linalg unchanged: trace of the 4-cycle adjacency is 0. */
    snprintf(buf, sizeof(buf), "Tr[AdjacencyMatrix[%s]]", dg);
    assert_eval_eq(buf, "0", 0);
    /* Det of the directed 4-cycle circulant is -1. */
    snprintf(buf, sizeof(buf), "Det[AdjacencyMatrix[%s]]", dg);
    assert_eval_eq(buf, "-1", 0);

    /* Round-trip: AdjacencyGraph[AdjacencyMatrix[g]] reproduces the edges. */
    snprintf(buf, sizeof(buf), "EdgeList[AdjacencyGraph[AdjacencyMatrix[%s]]]", dg);
    assert_eval_eq(buf, "{1 -> 2, 2 -> 3, 3 -> 4, 4 -> 1}", 0);
    /* Undirected round-trip stays undirected. */
    assert_eval_eq(
        "EdgeList[AdjacencyGraph[AdjacencyMatrix[Graph[{1,2,3},{1<->2,2<->3}]]]]",
        "{1 <-> 2, 2 <-> 3}", 0);

    /* Incidence matrix of an undirected path 1<->2<->3. */
    assert_eval_eq("IncidenceMatrix[Graph[{1,2,3},{1<->2,2<->3}]]",
                   "{{1, 0}, {1, 1}, {0, 1}}", 0);
}

/* ---- Phase 4: generators -------------------------------------------------- */
static void test_generators(void) {
    /* K5: 5 vertices, 10 edges, undirected, every vertex degree 4. */
    assert_eval_eq("VertexCount[CompleteGraph[5]]", "5", 0);
    assert_eval_eq("EdgeCount[CompleteGraph[5]]", "10", 0);
    assert_eval_eq("DirectedGraphQ[CompleteGraph[5]]", "False", 0);
    assert_eval_eq("VertexDegree[CompleteGraph[5]]", "{4, 4, 4, 4, 4}", 0);

    /* CycleGraph[n]: n vertices, n edges, every degree 2. */
    assert_eval_eq("EdgeCount[CycleGraph[5]]", "5", 0);
    assert_eval_eq("VertexDegree[CycleGraph[5]]", "{2, 2, 2, 2, 2}", 0);
    assert_eval_eq("EdgeList[CycleGraph[4]]",
                   "{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}", 0);

    /* PathGraph[n]: n vertices, n-1 edges; endpoints degree 1. */
    assert_eval_eq("EdgeCount[PathGraph[5]]", "4", 0);
    assert_eval_eq("VertexDegree[PathGraph[5]]", "{1, 2, 2, 2, 1}", 0);
    /* Explicit-vertex path. */
    assert_eval_eq("EdgeList[PathGraph[{a,b,c}]]", "{a <-> b, b <-> c}", 0);
}

static void test_random_graph(void) {
    /* Exactly m distinct edges, n vertices, valid simple graph. */
    assert_eval_eq("VertexCount[RandomGraph[{6, 5}]]", "6", 0);
    assert_eval_eq("EdgeCount[RandomGraph[{6, 5}]]", "5", 0);
    assert_eval_eq("GraphQ[RandomGraph[{6, 5}]]", "True", 0);
    /* Too many edges for a simple graph -> unevaluated (GraphQ False). */
    assert_eval_eq("GraphQ[RandomGraph[{3, 10}]]", "False", 0);
    /* Determinism under a fixed seed. */
    Expr* e1 = evaluate(parse_expression(
        "(SeedRandom[42]; EdgeList[RandomGraph[{6,5}]]) === "
        "(SeedRandom[42]; EdgeList[RandomGraph[{6,5}]])"));
    char* s = expr_to_string(e1);
    ASSERT(strcmp(s, "True") == 0);
    free(s);
    expr_free(e1);
}

/* ---- Phase 5: algorithms -------------------------------------------------- */
static void test_shortest_path(void) {
    /* Directed path 1->2->3->4. */
    const char* dg = "Graph[{1,2,3,4},{1->2,2->3,3->4}]";
    char buf[256];
    snprintf(buf, sizeof(buf), "FindShortestPath[%s, 1, 4]", dg);
    assert_eval_eq(buf, "{1, 2, 3, 4}", 0);
    snprintf(buf, sizeof(buf), "GraphDistance[%s, 1, 4]", dg);
    assert_eval_eq(buf, "3", 0);
    /* Direction matters: 4 cannot reach 1. */
    snprintf(buf, sizeof(buf), "GraphDistance[%s, 4, 1]", dg);
    assert_eval_eq(buf, "Infinity", 0);
    snprintf(buf, sizeof(buf), "FindShortestPath[%s, 4, 1]", dg);
    assert_eval_eq(buf, "{}", 0);
    /* Undirected: reachable both ways. */
    assert_eval_eq("GraphDistance[Graph[{1,2,3,4},{1<->2,2<->3,3<->4}], 4, 1]", "3", 0);
    /* A shortcut shortens the path. */
    assert_eval_eq(
        "GraphDistance[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4}], 1, 4]", "1", 0);
}

static void test_components(void) {
    /* Directed 4-cycle: one weak and one strong component. */
    assert_eval_eq("ConnectedComponents[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]",
                   "{{1, 2, 3, 4}}", 0);
    /* Two disjoint directed pieces -> two weak components. */
    assert_eval_eq("ConnectedComponents[Graph[{1,2,3,4},{1->2,3->4}]]",
                   "{{1, 2}, {3, 4}}", 0);
    /* Directed chain 1->2->3: weak = all together, strong = singletons. */
    assert_eval_eq("WeaklyConnectedComponents[Graph[{1,2,3},{1->2,2->3}]]",
                   "{{1, 2, 3}}", 0);
    assert_eval_eq("StronglyConnectedComponents[Graph[{1,2,3},{1->2,2->3}]]",
                   "{{1}, {2}, {3}}", 0);
    /* A directed cycle is one strong component. */
    assert_eval_eq("StronglyConnectedComponents[Graph[{1,2,3},{1->2,2->3,3->1}]]",
                   "{{1, 2, 3}}", 0);
}

static void test_spanning_and_connectivity(void) {
    /* Spanning tree of a connected graph has VertexCount - 1 edges. */
    assert_eval_eq("EdgeCount[FindSpanningTree[CompleteGraph[5]]]", "4", 0);
    assert_eval_eq("EdgeCount[FindSpanningTree[CycleGraph[6]]]", "5", 0);

    /* Connectivity. */
    assert_eval_eq("ConnectedGraphQ[PathGraph[4]]", "True", 0);
    assert_eval_eq("ConnectedGraphQ[Graph[{1,2,3,4},{1<->2,3<->4}]]", "False", 0);
    /* kappa: path=1, cycle=2, complete K4=3, disconnected=0. */
    assert_eval_eq("VertexConnectivity[PathGraph[4]]", "1", 0);
    assert_eval_eq("VertexConnectivity[CycleGraph[5]]", "2", 0);
    assert_eval_eq("VertexConnectivity[CompleteGraph[4]]", "3", 0);
    assert_eval_eq("VertexConnectivity[Graph[{1,2,3,4},{1<->2,3<->4}]]", "0", 0);
}

/* ---- Phase 6: GraphPlot --------------------------------------------------- */
static void test_graphplot(void) {
    /* Emits a Graphics object. */
    assert_eval_eq("Head[GraphPlot[CycleGraph[5]]]", "Graphics", 0);
    /* One Line per edge, one Disk per vertex. */
    assert_eval_eq("Count[GraphPlot[CycleGraph[5]], _Line, Infinity]", "5", 0);
    assert_eval_eq("Count[GraphPlot[CycleGraph[5]], _Disk, Infinity]", "5", 0);
    /* CompleteGraph[6]: 15 edges, 6 vertices. */
    assert_eval_eq("Count[GraphPlot[CompleteGraph[6]], _Line, Infinity]", "15", 0);
    assert_eval_eq("Count[GraphPlot[CompleteGraph[6]], _Disk, Infinity]", "6", 0);
    /* Non-graph argument stays unevaluated. */
    assert_eval_eq("Head[GraphPlot[5]]", "GraphPlot", 0);
}

static void test_graphplot_options(void) {
    /* Every layout still yields a Graphics with the same edge/vertex counts;
     * only the coordinates differ. */
    assert_eval_eq("Head[GraphPlot[CompleteGraph[6], GraphLayout -> \"SpringElectricalEmbedding\"]]", "Graphics", 0);
    assert_eval_eq("Count[GraphPlot[CompleteGraph[6], GraphLayout -> \"SpringElectricalEmbedding\"], _Line, Infinity]", "15", 0);
    assert_eval_eq("Count[GraphPlot[CycleGraph[8], GraphLayout -> \"GridEmbedding\"], _Disk, Infinity]", "8", 0);
    /* Unknown layout name falls back to circular (still valid). */
    assert_eval_eq("Head[GraphPlot[CycleGraph[5], GraphLayout -> \"NoSuchEmbedding\"]]", "Graphics", 0);
    /* VertexLabels -> None suppresses the Text labels; default draws them. */
    assert_eval_eq("Count[GraphPlot[CycleGraph[5], VertexLabels -> None], _Text, Infinity]", "0", 0);
    assert_eval_eq("Count[GraphPlot[CycleGraph[5]], _Text, Infinity]", "5", 0);
    /* Styling emits per-primitive RGBColor directives (one per edge + vertex). */
    assert_eval_eq("Count[GraphPlot[CycleGraph[5]], _RGBColor, Infinity]", "10", 0);
}

static void test_highlight_graph(void) {
    /* HighlightGraph returns a Graphics (not a Graph). */
    assert_eval_eq("Head[HighlightGraph[CycleGraph[8], {1, 2, 3}]]", "Graphics", 0);
    assert_eval_eq("Head[HighlightGraph[CompleteGraph[6], {1 <-> 2, 2 <-> 3}]]", "Graphics", 0);
    /* Edge/vertex counts are preserved (highlighting only recolors). */
    assert_eval_eq("Count[HighlightGraph[CycleGraph[6], {1, 2}], _Line, Infinity]", "6", 0);
    assert_eval_eq("Count[HighlightGraph[CycleGraph[6], {1, 2}], _Disk, Infinity]", "6", 0);
    /* A vertex list is treated as a path (accepted, still a Graphics). */
    assert_eval_eq("Head[HighlightGraph[PathGraph[5], {{1, 2, 3}}]]", "Graphics", 0);
    /* Non-graph / malformed arg stays unevaluated. */
    assert_eval_eq("Head[HighlightGraph[5, {1}]]", "HighlightGraph", 0);
}

static void test_graph3d(void) {
    /* Graph3D builds a canonical value with head Graph3D, and counts as a graph. */
    assert_eval_eq("GraphQ[Graph3D[{1,2},{1->2}]]", "True", 0);
    assert_eval_eq("Head[Graph3D[{1,2,3},{1<->2,2<->3}]]", "Graph3D", 0);
    /* Same vertex/edge readers work on a Graph3D value. */
    assert_eval_eq("VertexCount[Graph3D[{1,2,3},{1<->2,2<->3}]]", "3", 0);
    assert_eval_eq("EdgeCount[Graph3D[{1,2,3},{1<->2,2<->3}]]", "2", 0);
    /* Graph3D of an existing graph reuses its vertices/edges. */
    assert_eval_eq("VertexCount[Graph3D[CompleteGraph[5]]]", "5", 0);
    assert_eval_eq("EdgeCount[Graph3D[CompleteGraph[5]]]", "10", 0);
    /* InputForm round-trips the constructor. */
    assert_eval_eq("InputForm[Graph3D[{1,2},{1<->2}]]", "Graph3D[{1, 2}, {1 <-> 2}]", 0);
    /* Malformed (self-loop) stays unevaluated. */
    assert_eval_eq("Head[Graph3D[{1},{1->1}]]", "Graph3D", 0);
}

static void test_bipartite(void) {
    /* Even cycles / paths / stars / complete bipartite are 2-colorable. */
    assert_eval_eq("BipartiteGraphQ[CycleGraph[4]]", "True", 0);
    assert_eval_eq("BipartiteGraphQ[PathGraph[5]]", "True", 0);
    assert_eval_eq("BipartiteGraphQ[Graph[{0,1,2,3},{0<->1,0<->2,0<->3}]]", "True", 0);
    /* Odd cycles / triangles are not. */
    assert_eval_eq("BipartiteGraphQ[CycleGraph[5]]", "False", 0);
    assert_eval_eq("BipartiteGraphQ[CompleteGraph[3]]", "False", 0);
    /* Direction is ignored; edgeless is vacuously bipartite. */
    assert_eval_eq("BipartiteGraphQ[Graph[{1,2,3,4},{1->2,2->3,3->4,4->1}]]", "True", 0);
    assert_eval_eq("BipartiteGraphQ[Graph[{1,2,3},{}]]", "True", 0);
    /* One non-bipartite component makes the whole graph non-bipartite. */
    assert_eval_eq("BipartiteGraphQ[Graph[{1,2,3,4,5},{1<->2,2<->3,3<->1,4<->5}]]", "False", 0);
    /* Non-graph argument. */
    assert_eval_eq("BipartiteGraphQ[5]", "False", 0);
}

static void test_metrics(void) {
    /* Path P5: eccentricities {4,3,2,3,4}; diameter 4, radius 2, center {3}. */
    assert_eval_eq("VertexEccentricity[PathGraph[5]]", "{4, 3, 2, 3, 4}", 0);
    assert_eval_eq("GraphDiameter[PathGraph[5]]", "4", 0);
    assert_eval_eq("GraphRadius[PathGraph[5]]", "2", 0);
    assert_eval_eq("GraphCenter[PathGraph[5]]", "{3}", 0);
    /* Cycle C6: every eccentricity 3; center is all vertices. */
    assert_eval_eq("GraphDiameter[CycleGraph[6]]", "3", 0);
    assert_eval_eq("GraphCenter[CycleGraph[6]]", "{1, 2, 3, 4, 5, 6}", 0);
    /* Complete graph: diameter 1. */
    assert_eval_eq("GraphDiameter[CompleteGraph[5]]", "1", 0);
    /* Disconnected: diameter/radius Infinity, empty center. */
    assert_eval_eq("GraphDiameter[Graph[{1,2,3},{1<->2}]]", "Infinity", 0);
    assert_eval_eq("GraphRadius[Graph[{1,2,3},{1<->2}]]", "Infinity", 0);
    assert_eval_eq("GraphCenter[Graph[{1,2,3},{1<->2}]]", "{}", 0);
    /* Directed out-star: source reaches all (ecc 1), leaves unreachable →
       diameter Infinity but radius 1, center is the source. */
    assert_eval_eq("GraphRadius[Graph[{0,1,2,3},{0->1,0->2,0->3}]]", "1", 0);
    assert_eval_eq("GraphCenter[Graph[{0,1,2,3},{0->1,0->2,0->3}]]", "{0}", 0);
    assert_eval_eq("VertexEccentricity[Graph[{1},{}], 1]", "0", 0);
    /* Non-graph stays unevaluated. */
    assert_eval_eq("Head[GraphDiameter[5]]", "GraphDiameter", 0);
}

static void test_acyclic(void) {
    /* Directed DAG: a valid order + acyclic. */
    assert_eval_eq("TopologicalSort[Graph[{1,2,3},{1->2,2->3,1->3}]]", "{1, 2, 3}", 0);
    assert_eval_eq("AcyclicGraphQ[Graph[{1,2,3},{1->2,2->3,1->3}]]", "True", 0);
    /* Directed cycle: $Failed + not acyclic. */
    assert_eval_eq("TopologicalSort[Graph[{1,2,3},{1->2,2->3,3->1}]]", "$Failed", 0);
    assert_eval_eq("AcyclicGraphQ[Graph[{1,2,3},{1->2,2->3,3->1}]]", "False", 0);
    assert_eval_eq("AcyclicGraphQ[Graph[{1,2},{1->2,2->1}]]", "False", 0);
    /* Undirected forest is acyclic; an undirected cycle is not. Topological
       sort is $Failed for undirected edges (they act as 2-cycles). */
    assert_eval_eq("AcyclicGraphQ[PathGraph[4]]", "True", 0);
    assert_eval_eq("AcyclicGraphQ[Graph[{1,2,3,4,5},{1<->2,2<->3,4<->5}]]", "True", 0);
    assert_eval_eq("AcyclicGraphQ[CycleGraph[4]]", "False", 0);
    assert_eval_eq("TopologicalSort[PathGraph[4]]", "$Failed", 0);
    /* Edgeless graph: acyclic, order is all vertices. */
    assert_eval_eq("AcyclicGraphQ[Graph[{1,2,3},{}]]", "True", 0);
    assert_eval_eq("TopologicalSort[Graph[{1,2,3},{}]]", "{1, 2, 3}", 0);
    /* Non-graph. */
    assert_eval_eq("AcyclicGraphQ[5]", "False", 0);
}

static void test_complement(void) {
    /* Complement of a path keeps the one non-edge; complement of empty is K_n;
       complement of complete is edgeless. */
    assert_eval_eq("EdgeList[GraphComplement[PathGraph[3]]]", "{1 <-> 3}", 0);
    assert_eval_eq("EdgeCount[GraphComplement[Graph[{1,2,3},{}]]]", "3", 0);
    assert_eval_eq("EdgeCount[GraphComplement[CompleteGraph[4]]]", "0", 0);
    assert_eval_eq("EdgeCount[GraphComplement[CycleGraph[4]]]", "2", 0);
    /* Double complement restores the edge count. */
    assert_eval_eq("EdgeCount[GraphComplement[GraphComplement[CycleGraph[5]]]]", "5", 0);
    /* Directed complement stays directed (all ordered non-self pairs minus 1->2). */
    assert_eval_eq("EdgeCount[GraphComplement[Graph[{1,2,3},{1->2}]]]", "5", 0);
    assert_eval_eq("DirectedGraphQ[GraphComplement[Graph[{1,2,3},{1->2}]]]", "True", 0);
    assert_eval_eq("GraphQ[GraphComplement[PathGraph[4]]]", "True", 0);
}

static void test_kirchhoff(void) {
    /* Laplacian = D - A, exact small cases. */
    assert_eval_eq("KirchhoffMatrix[PathGraph[3]]", "{{1, -1, 0}, {-1, 2, -1}, {0, -1, 1}}", 0);
    assert_eval_eq("KirchhoffMatrix[CycleGraph[3]]", "{{2, -1, -1}, {-1, 2, -1}, {-1, -1, 2}}", 0);
    /* Every row sums to 0; trace = sum of degrees = 2|E|. */
    assert_eval_eq("Total[KirchhoffMatrix[CycleGraph[4]]]", "{0, 0, 0, 0}", 0);
    assert_eval_eq("Tr[KirchhoffMatrix[CycleGraph[4]]]", "8", 0);
    /* Interop: Laplacian eigenvalues include 0 (once, since connected). */
    assert_eval_eq("Eigenvalues[KirchhoffMatrix[PathGraph[3]]]", "{3, 1, 0}", 0);
    assert_eval_eq("Head[KirchhoffMatrix[5]]", "KirchhoffMatrix", 0);
}

static void test_edge_connectivity(void) {
    /* Complete K_n -> n-1; cycle -> 2; path/tree/star -> 1. */
    assert_eval_eq("EdgeConnectivity[CompleteGraph[4]]", "3", 0);
    assert_eval_eq("EdgeConnectivity[CycleGraph[5]]", "2", 0);
    assert_eval_eq("EdgeConnectivity[PathGraph[4]]", "1", 0);
    assert_eval_eq("EdgeConnectivity[Graph[{0,1,2,3},{0<->1,0<->2,0<->3}]]", "1", 0);
    /* A bridge between two triangles -> 1. */
    assert_eval_eq("EdgeConnectivity[Graph[{1,2,3,4,5,6},{1<->2,2<->3,3<->1,4<->5,5<->6,6<->4,3<->4}]]", "1", 0);
    /* Disconnected -> 0. */
    assert_eval_eq("EdgeConnectivity[Graph[{1,2,3},{1<->2}]]", "0", 0);
    /* Directed: strongly connected cycle -> 1; a path is not, -> 0. */
    assert_eval_eq("EdgeConnectivity[Graph[{1,2,3},{1->2,2->3,3->1}]]", "1", 0);
    assert_eval_eq("EdgeConnectivity[Graph[{1,2,3},{1->2,2->3}]]", "0", 0);
    assert_eval_eq("Head[EdgeConnectivity[5]]", "EdgeConnectivity", 0);
}

static void test_line_graph(void) {
    /* L(P3): the two edges share vertex 2 → one edge; vertices ARE the edges. */
    assert_eval_eq("VertexCount[LineGraph[PathGraph[3]]]", "2", 0);
    assert_eval_eq("EdgeCount[LineGraph[PathGraph[3]]]", "1", 0);
    assert_eval_eq("VertexList[LineGraph[PathGraph[3]]]", "{1 <-> 2, 2 <-> 3}", 0);
    /* L(C3) = K3; L(Cn) = Cn; L(K4) has 6 vertices and 12 edges. */
    assert_eval_eq("EdgeCount[LineGraph[CycleGraph[3]]]", "3", 0);
    assert_eval_eq("EdgeCount[LineGraph[CycleGraph[5]]]", "5", 0);
    assert_eval_eq("VertexCount[LineGraph[CompleteGraph[4]]]", "6", 0);
    assert_eval_eq("EdgeCount[LineGraph[CompleteGraph[4]]]", "12", 0);
    /* Star K1,3 → the three edges pairwise share the center → K3. */
    assert_eval_eq("EdgeCount[LineGraph[Graph[{0,1,2,3},{0<->1,0<->2,0<->3}]]]", "3", 0);
    /* Directed line graph joins arcs head-to-tail and stays directed. */
    assert_eval_eq("EdgeCount[LineGraph[Graph[{1,2,3},{1->2,2->3}]]]", "1", 0);
    assert_eval_eq("DirectedGraphQ[LineGraph[Graph[{1,2,3},{1->2,2->3}]]]", "True", 0);
}

static void test_eulerian(void) {
    /* Even + connected → Eulerian; odd degree → not. */
    assert_eval_eq("EulerianGraphQ[CycleGraph[5]]", "True", 0);
    assert_eval_eq("EulerianGraphQ[PathGraph[4]]", "False", 0);
    assert_eval_eq("EulerianGraphQ[CompleteGraph[3]]", "True", 0);
    assert_eval_eq("EulerianGraphQ[CompleteGraph[4]]", "False", 0);
    assert_eval_eq("EulerianGraphQ[CompleteGraph[5]]", "True", 0);
    /* Figure-eight (two triangles sharing a vertex) is Eulerian; two disjoint
       cycles are not (disconnected). An isolated vertex is ignored. */
    assert_eval_eq("EulerianGraphQ[Graph[{1,2,3,4,5},{1<->2,2<->3,3<->1,1<->4,4<->5,5<->1}]]", "True", 0);
    assert_eval_eq("EulerianGraphQ[Graph[{1,2,3,4,5,6},{1<->2,2<->3,3<->1,4<->5,5<->6,6<->4}]]", "False", 0);
    assert_eval_eq("EulerianGraphQ[Graph[{1,2,3,4},{1<->2,2<->3,3<->1}]]", "True", 0);
    assert_eval_eq("EulerianGraphQ[Graph[{1,2,3},{}]]", "True", 0);
    /* Directed: in==out everywhere and connected → Eulerian; a path is not. */
    assert_eval_eq("EulerianGraphQ[Graph[{1,2,3},{1->2,2->3,3->1}]]", "True", 0);
    assert_eval_eq("EulerianGraphQ[Graph[{1,2,3},{1->2,2->3}]]", "False", 0);
    assert_eval_eq("EulerianGraphQ[5]", "False", 0);
}

static void test_star_wheel(void) {
    /* Star K_{1,n-1}: center joined to n-1 leaves; bipartite. */
    assert_eval_eq("VertexCount[StarGraph[5]]", "5", 0);
    assert_eval_eq("EdgeCount[StarGraph[5]]", "4", 0);
    assert_eval_eq("VertexDegree[StarGraph[5], 1]", "4", 0);
    assert_eval_eq("BipartiteGraphQ[StarGraph[5]]", "True", 0);
    assert_eval_eq("EdgeList[StarGraph[4]]", "{1 <-> 2, 1 <-> 3, 1 <-> 4}", 0);
    /* Wheel: rim cycle + hub; 2(n-1) edges, hub degree n-1, has triangles so
       it is not bipartite; W_4 = K_4. Small n stays unevaluated. */
    assert_eval_eq("VertexCount[WheelGraph[5]]", "5", 0);
    assert_eval_eq("EdgeCount[WheelGraph[5]]", "8", 0);
    assert_eval_eq("VertexDegree[WheelGraph[5], 5]", "4", 0);
    assert_eval_eq("BipartiteGraphQ[WheelGraph[5]]", "False", 0);
    assert_eval_eq("EdgeCount[WheelGraph[4]]", "6", 0);
    assert_eval_eq("Head[WheelGraph[3]]", "WheelGraph", 0);
}

static void test_complete_multipartite(void) {
    /* CompleteGraph[{m,n}] is complete bipartite K_{m,n}. */
    assert_eval_eq("VertexCount[CompleteGraph[{2,3}]]", "5", 0);
    assert_eval_eq("EdgeCount[CompleteGraph[{2,3}]]", "6", 0);
    assert_eval_eq("BipartiteGraphQ[CompleteGraph[{2,3}]]", "True", 0);
    assert_eval_eq("EdgeCount[CompleteGraph[{3,3}]]", "9", 0);
    assert_eval_eq("VertexDegree[CompleteGraph[{2,3}], 1]", "3", 0);
    /* Complete multipartite: K_{2,2,2} is the octahedron (12 edges, not bip.). */
    assert_eval_eq("EdgeCount[CompleteGraph[{2,2,2}]]", "12", 0);
    assert_eval_eq("BipartiteGraphQ[CompleteGraph[{2,2,2}]]", "False", 0);
    /* The integer form is unchanged. */
    assert_eval_eq("EdgeCount[CompleteGraph[5]]", "10", 0);
}

static void test_grid_hypercube(void) {
    /* Grid: 2x3 has 7 edges, k-dim works, {n} is a path, and it is bipartite. */
    assert_eval_eq("VertexCount[GridGraph[{2,3}]]", "6", 0);
    assert_eval_eq("EdgeCount[GridGraph[{2,3}]]", "7", 0);
    assert_eval_eq("EdgeCount[GridGraph[{3,3}]]", "12", 0);
    assert_eval_eq("BipartiteGraphQ[GridGraph[{3,3}]]", "True", 0);
    assert_eval_eq("EdgeCount[GridGraph[{4}]]", "3", 0);
    assert_eval_eq("Head[GridGraph[5]]", "GridGraph", 0);
    /* Hypercube Q_k: 2^k vertices, k*2^(k-1) edges, k-regular, bipartite;
       Q_2 = C_4; the 2x2x2 grid is Q_3. */
    assert_eval_eq("VertexCount[HypercubeGraph[3]]", "8", 0);
    assert_eval_eq("EdgeCount[HypercubeGraph[3]]", "12", 0);
    assert_eval_eq("BipartiteGraphQ[HypercubeGraph[3]]", "True", 0);
    assert_eval_eq("VertexDegree[HypercubeGraph[3], 1]", "3", 0);
    assert_eval_eq("EdgeCount[HypercubeGraph[2]]", "4", 0);
    assert_eval_eq("VertexCount[HypercubeGraph[0]]", "1", 0);
    assert_eval_eq("EdgeCount[GridGraph[{2,2,2}]]", "12", 0);
}

static void test_closeness(void) {
    /* Exact rationals: complete → all 1; cycle C4 → all 3/4; path P3 → ends 2/3,
       middle 1; star center 1, leaves 3/5. */
    assert_eval_eq("ClosenessCentrality[CompleteGraph[4]]", "{1, 1, 1, 1}", 0);
    assert_eval_eq("ClosenessCentrality[CycleGraph[4]]", "{3/4, 3/4, 3/4, 3/4}", 0);
    assert_eval_eq("ClosenessCentrality[PathGraph[3]]", "{2/3, 1, 2/3}", 0);
    assert_eval_eq("ClosenessCentrality[StarGraph[4]]", "{1, 3/5, 3/5, 3/5}", 0);
    /* Isolated / disconnected vertices score 0. */
    assert_eval_eq("ClosenessCentrality[Graph[{1,2},{}]]", "{0, 0}", 0);
    /* Directed: distances follow direction; a sink reaches nobody → 0. */
    assert_eval_eq("ClosenessCentrality[Graph[{1,2,3},{1->2,2->3}]]", "{2/3, 1/2, 0}", 0);
    assert_eval_eq("Head[ClosenessCentrality[5]]", "ClosenessCentrality", 0);
}

static void test_transitive_closure(void) {
    /* Directed chain: closure adds the reachable 1->3; stays directed. */
    assert_eval_eq("EdgeList[TransitiveClosure[Graph[{1,2,3},{1->2,2->3}]]]", "{1 -> 2, 1 -> 3, 2 -> 3}", 0);
    assert_eval_eq("DirectedGraphQ[TransitiveClosure[Graph[{1,2,3},{1->2,2->3}]]]", "True", 0);
    /* Directed cycle: everyone reaches everyone → complete digraph (6 arcs). */
    assert_eval_eq("EdgeCount[TransitiveClosure[Graph[{1,2,3},{1->2,2->3,3->1}]]]", "6", 0);
    /* Undirected: each component becomes a clique. */
    assert_eval_eq("EdgeList[TransitiveClosure[PathGraph[3]]]", "{1 <-> 2, 1 <-> 3, 2 <-> 3}", 0);
    assert_eval_eq("EdgeCount[TransitiveClosure[CycleGraph[4]]]", "6", 0);
    assert_eval_eq("EdgeCount[TransitiveClosure[Graph[{1,2,3,4,5},{1<->2,3<->4,4<->5}]]]", "4", 0);
    assert_eval_eq("EdgeCount[TransitiveClosure[CompleteGraph[4]]]", "6", 0);
    assert_eval_eq("Head[TransitiveClosure[5]]", "TransitiveClosure", 0);
}

static void test_betweenness(void) {
    /* Path: internal vertex k has (k-1)(n-k); complete has all 0; star center
       carries every leaf-leaf path. */
    assert_eval_eq("BetweennessCentrality[PathGraph[3]]", "{0, 1, 0}", 0);
    assert_eval_eq("BetweennessCentrality[PathGraph[5]]", "{0, 3, 4, 3, 0}", 0);
    assert_eval_eq("BetweennessCentrality[CompleteGraph[4]]", "{0, 0, 0, 0}", 0);
    assert_eval_eq("BetweennessCentrality[StarGraph[5]]", "{6, 0, 0, 0, 0}", 0);
    /* Tied shortest paths give exact fractions: C4 antipodal pairs → 1/2 each. */
    assert_eval_eq("BetweennessCentrality[CycleGraph[4]]", "{1/2, 1/2, 1/2, 1/2}", 0);
    /* Directed counts ordered pairs (no halving). */
    assert_eval_eq("BetweennessCentrality[Graph[{1,2,3},{1->2,2->3}]]", "{0, 1, 0}", 0);
    assert_eval_eq("BetweennessCentrality[Graph[{1,2,3},{}]]", "{0, 0, 0}", 0);
    assert_eval_eq("Head[BetweennessCentrality[5]]", "BetweennessCentrality", 0);
}

static void test_find_eulerian(void) {
    /* Eulerian graphs return a closed tour of ne+1 vertices. */
    assert_eval_eq("FindEulerianCycle[CycleGraph[3]]", "{1, 2, 3, 1}", 0);
    assert_eval_eq("Length[FindEulerianCycle[CycleGraph[4]]]", "5", 0);
    assert_eval_eq("Length[FindEulerianCycle[CompleteGraph[5]]]", "11", 0);
    assert_eval_eq("First[FindEulerianCycle[CycleGraph[4]]] === Last[FindEulerianCycle[CycleGraph[4]]]", "True", 0);
    /* Directed Eulerian circuit. */
    assert_eval_eq("Length[FindEulerianCycle[Graph[{1,2,3},{1->2,2->3,3->1}]]]", "4", 0);
    /* Non-Eulerian (odd degree, disconnected, open directed, edgeless) → {}. */
    assert_eval_eq("FindEulerianCycle[PathGraph[3]]", "{}", 0);
    assert_eval_eq("FindEulerianCycle[CompleteGraph[4]]", "{}", 0);
    assert_eval_eq("FindEulerianCycle[Graph[{1,2,3,4,5,6},{1<->2,2<->3,3<->1,4<->5,5<->6,6<->4}]]", "{}", 0);
    assert_eval_eq("FindEulerianCycle[Graph[{1,2,3},{1->2,2->3}]]", "{}", 0);
    assert_eval_eq("FindEulerianCycle[Graph[{1,2},{}]]", "{}", 0);
    assert_eval_eq("Head[FindEulerianCycle[5]]", "FindEulerianCycle", 0);
}

static void test_find_hamiltonian(void) {
    /* Cycles/complete graphs have Hamiltonian cycles: n+1 closed vertex list. */
    assert_eval_eq("FindHamiltonianCycle[CycleGraph[3]]", "{1, 2, 3, 1}", 0);
    assert_eval_eq("Length[FindHamiltonianCycle[CycleGraph[5]]]", "6", 0);
    assert_eval_eq("Length[FindHamiltonianCycle[CompleteGraph[4]]]", "5", 0);
    assert_eval_eq("Sort[Union[FindHamiltonianCycle[CompleteGraph[4]]]]", "{1, 2, 3, 4}", 0);
    assert_eval_eq("First[FindHamiltonianCycle[CycleGraph[5]]] === Last[FindHamiltonianCycle[CycleGraph[5]]]", "True", 0);
    assert_eval_eq("Length[FindHamiltonianCycle[WheelGraph[5]]]", "6", 0);
    /* Directed Hamiltonian cycle follows arc direction. */
    assert_eval_eq("Length[FindHamiltonianCycle[Graph[{1,2,3},{1->2,2->3,3->1}]]]", "4", 0);
    /* No Hamiltonian cycle: paths, stars, n<3, edgeless, disconnected, broken directed → {}. */
    assert_eval_eq("FindHamiltonianCycle[PathGraph[4]]", "{}", 0);
    assert_eval_eq("FindHamiltonianCycle[StarGraph[5]]", "{}", 0);
    assert_eval_eq("FindHamiltonianCycle[Graph[{1,2},{1<->2}]]", "{}", 0);
    assert_eval_eq("FindHamiltonianCycle[Graph[{1,2,3},{}]]", "{}", 0);
    assert_eval_eq("FindHamiltonianCycle[Graph[{1,2,3,4,5,6},{1<->2,2<->3,3<->1,4<->5,5<->6,6<->4}]]", "{}", 0);
    assert_eval_eq("FindHamiltonianCycle[Graph[{1,2,3},{1->2,2->3,1->3}]]", "{}", 0);
    assert_eval_eq("Head[FindHamiltonianCycle[7]]", "FindHamiltonianCycle", 0);
}

static void test_graph_power(void) {
    /* Undirected: distance-<=k closure. P4^2 adds {1,3},{2,4}; P4^3 = K4. */
    assert_eval_eq("EdgeCount[GraphPower[PathGraph[4], 2]]", "5", 0);
    assert_eval_eq("EdgeCount[GraphPower[PathGraph[4], 3]]", "6", 0);
    assert_eval_eq("EdgeCount[GraphPower[CycleGraph[5], 2]]", "10", 0);  /* K5 */
    assert_eval_eq("EdgeCount[GraphPower[CycleGraph[5], 1]]", "5", 0);   /* unchanged */
    assert_eval_eq("VertexList[GraphPower[PathGraph[4], 2]]", "{1, 2, 3, 4}", 0);
    assert_eval_eq("FreeQ[EdgeList[GraphPower[CycleGraph[4], 2]], UndirectedEdge[x_, x_]]", "True", 0);
    assert_eval_eq("EdgeCount[GraphPower[Graph[{1,2,3},{}], 2]]", "0", 0);
    /* Directed: reachability follows arcs; power stays directed. */
    assert_eval_eq("DirectedGraphQ[GraphPower[Graph[{1,2,3},{1->2,2->3}], 2]]", "True", 0);
    assert_eval_eq("MemberQ[EdgeList[GraphPower[Graph[{1,2,3},{1->2,2->3}], 2]], DirectedEdge[1,3]]", "True", 0);
    /* Non-positive / symbolic k stays unevaluated. */
    assert_eval_eq("Head[GraphPower[CycleGraph[3], 0]]", "GraphPower", 0);
    assert_eval_eq("Head[GraphPower[CycleGraph[3], k]]", "GraphPower", 0);
}

static void test_find_cycle(void) {
    /* Undirected: a cycle as {{edges...}}, or {} when acyclic. */
    assert_eval_eq("FindCycle[CycleGraph[3]]", "{{1 <-> 2, 2 <-> 3, 3 <-> 1}}", 0);
    assert_eval_eq("Length[First[FindCycle[CycleGraph[4]]]]", "4", 0);
    assert_eval_eq("MatchQ[First[FindCycle[CycleGraph[4]]], {UndirectedEdge[_,_]..}]", "True", 0);
    assert_eval_eq("FindCycle[CompleteGraph[5]] =!= {}", "True", 0);
    assert_eval_eq("FindCycle[PathGraph[5]]", "{}", 0);
    assert_eval_eq("FindCycle[Graph[{1,2,3,4},{1<->2,2<->3,3<->4}]]", "{}", 0);
    assert_eval_eq("FindCycle[Graph[{1,2,3},{}]]", "{}", 0);
    /* The returned edges close up into a cycle. */
    assert_eval_eq("With[{c=First[FindCycle[CompleteGraph[5]]]}, First[First[c]] === Last[Last[c]]]", "True", 0);
    /* Directed: cycles follow arcs; DAGs give {}. */
    assert_eval_eq("FindCycle[Graph[{1,2,3},{1->2,2->3,3->1}]]", "{{1 -> 2, 2 -> 3, 3 -> 1}}", 0);
    assert_eval_eq("FindCycle[Graph[{1,2},{1->2,2->1}]]", "{{1 -> 2, 2 -> 1}}", 0);
    assert_eval_eq("FindCycle[Graph[{1,2,3},{1->2,2->3,1->3}]]", "{}", 0);
    assert_eval_eq("Head[FindCycle[9]]", "FindCycle", 0);
}

static void test_graph_distance_matrix(void) {
    /* Undirected: symmetric, zero diagonal, BFS distances. */
    assert_eval_eq("GraphDistanceMatrix[PathGraph[3]]", "{{0, 1, 2}, {1, 0, 1}, {2, 1, 0}}", 0);
    assert_eval_eq("GraphDistanceMatrix[CompleteGraph[3]]", "{{0, 1, 1}, {1, 0, 1}, {1, 1, 0}}", 0);
    assert_eval_eq("Tr[GraphDistanceMatrix[CompleteGraph[5]]]", "0", 0);
    assert_eval_eq("With[{m=GraphDistanceMatrix[CycleGraph[5]]}, m === Transpose[m]]", "True", 0);
    assert_eval_eq("GraphDistanceMatrix[CycleGraph[6]][[1,4]]", "3", 0);
    /* Unreachable pairs are Infinity. */
    assert_eval_eq("GraphDistanceMatrix[Graph[{1,2,3},{1<->2}]]", "{{0, 1, Infinity}, {1, 0, Infinity}, {Infinity, Infinity, 0}}", 0);
    assert_eval_eq("GraphDistanceMatrix[Graph[{1,2},{}]]", "{{0, Infinity}, {Infinity, 0}}", 0);
    assert_eval_eq("GraphDistanceMatrix[Graph[{1},{}]]", "{{0}}", 0);
    /* Directed: follows arc direction (asymmetric). */
    assert_eval_eq("GraphDistanceMatrix[Graph[{1,2,3},{1->2,2->3}]]", "{{0, 1, 2}, {Infinity, 0, 1}, {Infinity, Infinity, 0}}", 0);
    /* Agrees with the pairwise GraphDistance builtin. */
    assert_eval_eq("GraphDistanceMatrix[PathGraph[4]][[1,4]] == GraphDistance[PathGraph[4],1,4]", "True", 0);
    assert_eval_eq("Head[GraphDistanceMatrix[5]]", "GraphDistanceMatrix", 0);
}

static void test_graph_density(void) {
    /* Undirected: 2m / (n(n-1)), reduced exact rational. */
    assert_eval_eq("GraphDensity[CompleteGraph[5]]", "1", 0);
    assert_eval_eq("GraphDensity[Graph[{1,2,3,4},{}]]", "0", 0);
    assert_eval_eq("GraphDensity[CycleGraph[4]]", "2/3", 0);
    assert_eval_eq("GraphDensity[PathGraph[4]]", "1/2", 0);
    assert_eval_eq("GraphDensity[StarGraph[5]]", "2/5", 0);
    assert_eval_eq("GraphDensity[Graph[{1,2},{1<->2}]]", "1", 0);
    assert_eval_eq("GraphDensity[Graph[{1},{}]]", "0", 0);
    /* Directed: m / (n(n-1)); a complete digraph is 1. */
    assert_eval_eq("GraphDensity[Graph[{1,2,3},{1->2,2->3}]]", "1/3", 0);
    assert_eval_eq("GraphDensity[Graph[{1,2,3},{1->2,2->1,1->3,3->1,2->3,3->2}]]", "1", 0);
    assert_eval_eq("Head[GraphDensity[7]]", "GraphDensity", 0);
}

static void test_degree_centrality(void) {
    assert_eval_eq("DegreeCentrality[CycleGraph[4]]", "{2, 2, 2, 2}", 0);
    assert_eval_eq("DegreeCentrality[CompleteGraph[4]]", "{3, 3, 3, 3}", 0);
    assert_eval_eq("DegreeCentrality[StarGraph[5]]", "{4, 1, 1, 1, 1}", 0);
    assert_eval_eq("DegreeCentrality[PathGraph[4]]", "{1, 2, 2, 1}", 0);
    assert_eval_eq("DegreeCentrality[Graph[{1,2,3},{}]]", "{0, 0, 0}", 0);
    /* Directed: in-degree + out-degree. */
    assert_eval_eq("DegreeCentrality[Graph[{1,2,3},{1->2,2->3}]]", "{1, 2, 1}", 0);
    /* Consistency: matches VertexDegree (undirected) and obeys the handshake lemma. */
    assert_eval_eq("DegreeCentrality[CycleGraph[5]] == VertexDegree[CycleGraph[5]]", "True", 0);
    assert_eval_eq("Total[DegreeCentrality[CompleteGraph[6]]] == 2*EdgeCount[CompleteGraph[6]]", "True", 0);
    assert_eval_eq("Head[DegreeCentrality[5]]", "DegreeCentrality", 0);
}

static void test_find_hamiltonian_path(void) {
    assert_eval_eq("FindHamiltonianPath[PathGraph[4]]", "{1, 2, 3, 4}", 0);
    assert_eval_eq("Length[FindHamiltonianPath[CycleGraph[4]]]", "4", 0);
    assert_eval_eq("Sort[FindHamiltonianPath[CompleteGraph[4]]]", "{1, 2, 3, 4}", 0);
    assert_eval_eq("With[{p=FindHamiltonianPath[CompleteGraph[5]]}, Length[Union[p]] == 5]", "True", 0);
    assert_eval_eq("FindHamiltonianPath[Graph[{1},{}]]", "{1}", 0);
    /* No Hamiltonian path: big star, edgeless, disconnected → {}. */
    assert_eval_eq("FindHamiltonianPath[StarGraph[5]]", "{}", 0);
    assert_eval_eq("FindHamiltonianPath[Graph[{1,2,3},{}]]", "{}", 0);
    assert_eval_eq("FindHamiltonianPath[Graph[{1,2,3,4},{1<->2,3<->4}]]", "{}", 0);
    /* Directed follows arcs. */
    assert_eval_eq("FindHamiltonianPath[Graph[{1,2,3},{1->2,2->3}]]", "{1, 2, 3}", 0);
    assert_eval_eq("FindHamiltonianPath[Graph[{1,2,3},{1->2,3->2}]]", "{}", 0);
    /* A path graph has a Hamiltonian path but no Hamiltonian cycle. */
    assert_eval_eq("FindHamiltonianPath[PathGraph[4]] =!= {} && FindHamiltonianCycle[PathGraph[4]] === {}", "True", 0);
    assert_eval_eq("Head[FindHamiltonianPath[5]]", "FindHamiltonianPath", 0);
}

static void test_kcore_components(void) {
    /* Complete/cycle cores: whole graph at the degree, empty just above. */
    assert_eval_eq("KCoreComponents[CompleteGraph[4], 3]", "{{1, 2, 3, 4}}", 0);
    assert_eval_eq("KCoreComponents[CompleteGraph[4], 4]", "{}", 0);
    assert_eval_eq("KCoreComponents[CycleGraph[5], 2]", "{{1, 2, 3, 4, 5}}", 0);
    assert_eval_eq("KCoreComponents[CycleGraph[5], 3]", "{}", 0);
    /* Cascade peeling: a path collapses entirely at k = 2. */
    assert_eval_eq("KCoreComponents[PathGraph[4], 1]", "{{1, 2, 3, 4}}", 0);
    assert_eval_eq("KCoreComponents[PathGraph[4], 2]", "{}", 0);
    /* k = 0 keeps every vertex; edgeless → singletons, empty at k = 1. */
    assert_eval_eq("KCoreComponents[Graph[{1,2,3},{}], 0]", "{{1}, {2}, {3}}", 0);
    assert_eval_eq("KCoreComponents[Graph[{1,2,3},{}], 1]", "{}", 0);
    /* Multiple components and pendant removal. */
    assert_eval_eq("KCoreComponents[Graph[{1,2,3,4,5,6},{1<->2,2<->3,3<->1,4<->5,5<->6,6<->4}], 2]", "{{1, 2, 3}, {4, 5, 6}}", 0);
    assert_eval_eq("KCoreComponents[Graph[{1,2,3,4},{1<->2,2<->3,3<->1,3<->4}], 2]", "{{1, 2, 3}}", 0);
    /* Direction ignored (underlying undirected graph). */
    assert_eval_eq("KCoreComponents[Graph[{1,2,3},{1->2,2->3,3->1}], 2]", "{{1, 2, 3}}", 0);
    assert_eval_eq("Head[KCoreComponents[CycleGraph[3], k]]", "KCoreComponents", 0);
    assert_eval_eq("Head[KCoreComponents[5, 2]]", "KCoreComponents", 0);
}

static void test_local_clustering(void) {
    /* Cliques → 1 everywhere; triangle-free (cycles, stars, paths) → 0. */
    assert_eval_eq("LocalClusteringCoefficient[CompleteGraph[4]]", "{1, 1, 1, 1}", 0);
    assert_eval_eq("LocalClusteringCoefficient[CompleteGraph[3]]", "{1, 1, 1}", 0);
    assert_eval_eq("LocalClusteringCoefficient[CycleGraph[5]]", "{0, 0, 0, 0, 0}", 0);
    assert_eval_eq("LocalClusteringCoefficient[PathGraph[3]]", "{0, 0, 0}", 0);
    assert_eval_eq("LocalClusteringCoefficient[StarGraph[5]]", "{0, 0, 0, 0, 0}", 0);
    assert_eval_eq("LocalClusteringCoefficient[Graph[{1,2,3},{}]]", "{0, 0, 0}", 0);
    /* Fractional cases: exact rationals. */
    assert_eval_eq("LocalClusteringCoefficient[Graph[{1,2,3,4},{1<->2,2<->3,3<->1,3<->4}]]", "{1, 1, 1/3, 0}", 0);
    assert_eval_eq("LocalClusteringCoefficient[Graph[{1,2,3,4},{1<->2,1<->3,1<->4,2<->3,2<->4}]]", "{2/3, 2/3, 1, 1}", 0);
    /* Direction ignored (underlying undirected graph). */
    assert_eval_eq("LocalClusteringCoefficient[Graph[{1,2,3},{1->2,2->3,3->1}]]", "{1, 1, 1}", 0);
    assert_eval_eq("Head[LocalClusteringCoefficient[5]]", "LocalClusteringCoefficient", 0);
}

static void test_global_clustering(void) {
    /* Cliques → 1; triangle-free → 0. */
    assert_eval_eq("GlobalClusteringCoefficient[CompleteGraph[4]]", "1", 0);
    assert_eval_eq("GlobalClusteringCoefficient[CompleteGraph[3]]", "1", 0);
    assert_eval_eq("GlobalClusteringCoefficient[CycleGraph[5]]", "0", 0);
    assert_eval_eq("GlobalClusteringCoefficient[PathGraph[3]]", "0", 0);
    assert_eval_eq("GlobalClusteringCoefficient[StarGraph[6]]", "0", 0);
    assert_eval_eq("GlobalClusteringCoefficient[Graph[{1,2},{1<->2}]]", "0", 0);
    assert_eval_eq("GlobalClusteringCoefficient[Graph[{1,2,3},{}]]", "0", 0);
    /* Fractional transitivity: exact rationals. */
    assert_eval_eq("GlobalClusteringCoefficient[Graph[{1,2,3,4},{1<->2,2<->3,3<->1,3<->4}]]", "3/5", 0);
    assert_eval_eq("GlobalClusteringCoefficient[Graph[{1,2,3,4},{1<->2,1<->3,1<->4,2<->3,2<->4}]]", "3/4", 0);
    /* Direction ignored. */
    assert_eval_eq("GlobalClusteringCoefficient[Graph[{1,2,3},{1->2,2->3,3->1}]]", "1", 0);
    assert_eval_eq("Head[GlobalClusteringCoefficient[5]]", "GlobalClusteringCoefficient", 0);
}

static void test_mean_clustering(void) {
    assert_eval_eq("MeanClusteringCoefficient[CompleteGraph[4]]", "1", 0);
    assert_eval_eq("MeanClusteringCoefficient[CycleGraph[5]]", "0", 0);
    assert_eval_eq("MeanClusteringCoefficient[StarGraph[6]]", "0", 0);
    assert_eval_eq("MeanClusteringCoefficient[Graph[{1,2,3},{}]]", "0", 0);
    assert_eval_eq("MeanClusteringCoefficient[Graph[{1,2,3,4},{1<->2,2<->3,3<->1,3<->4}]]", "7/12", 0);
    assert_eval_eq("MeanClusteringCoefficient[Graph[{1,2,3,4},{1<->2,1<->3,1<->4,2<->3,2<->4}]]", "5/6", 0);
    /* Equals the mean of the local coefficients; differs from global transitivity. */
    assert_eval_eq("MeanClusteringCoefficient[Graph[{1,2,3,4},{1<->2,2<->3,3<->1,3<->4}]] == Total[LocalClusteringCoefficient[Graph[{1,2,3,4},{1<->2,2<->3,3<->1,3<->4}]]]/4", "True", 0);
    assert_eval_eq("MeanClusteringCoefficient[Graph[{1,2,3,4},{1<->2,1<->3,1<->4,2<->3,2<->4}]] != GlobalClusteringCoefficient[Graph[{1,2,3,4},{1<->2,1<->3,1<->4,2<->3,2<->4}]]", "True", 0);
    assert_eval_eq("MeanClusteringCoefficient[Graph[{1,2,3},{1->2,2->3,3->1}]]", "1", 0);
    assert_eval_eq("Head[MeanClusteringCoefficient[5]]", "MeanClusteringCoefficient", 0);
}

static void test_find_clique(void) {
    /* Cliques of complete graphs are the whole graph. */
    assert_eval_eq("FindClique[CompleteGraph[4]]", "{{1, 2, 3, 4}}", 0);
    assert_eval_eq("FindClique[CompleteGraph[3]]", "{{1, 2, 3}}", 0);
    assert_eval_eq("FindClique[Graph[{1,2,3},{1<->2,2<->3,3<->1}]]", "{{1, 2, 3}}", 0);
    assert_eval_eq("FindClique[Graph[{1,2},{1<->2}]]", "{{1, 2}}", 0);
    /* Max clique sizes. */
    assert_eval_eq("Length[First[FindClique[CycleGraph[5]]]]", "2", 0);
    assert_eval_eq("Length[First[FindClique[CycleGraph[4]]]]", "2", 0);
    assert_eval_eq("Length[First[FindClique[Graph[{1,2,3,4},{1<->2,1<->3,1<->4,2<->3,2<->4}]]]]", "3", 0);
    assert_eval_eq("Length[First[FindClique[Graph[{1,2,3,4,5,6},{1<->2,2<->3,3<->1,4<->5,5<->6,6<->4}]]]]", "3", 0);
    assert_eval_eq("Length[First[FindClique[Graph[{1,2,3},{}]]]]", "1", 0);
    /* Direction ignored. */
    assert_eval_eq("Length[First[FindClique[Graph[{1,2,3},{1->2,2->3,3->1}]]]]", "3", 0);
    assert_eval_eq("Head[FindClique[5]]", "FindClique", 0);
}

static void test_find_independent(void) {
    /* Edgeless → all vertices; complete → a singleton. */
    assert_eval_eq("FindIndependentVertexSet[Graph[{1,2,3},{}]]", "{{1, 2, 3}}", 0);
    assert_eval_eq("Length[First[FindIndependentVertexSet[CompleteGraph[4]]]]", "1", 0);
    /* Cycles / paths / star sizes. */
    assert_eval_eq("FindIndependentVertexSet[CycleGraph[4]]", "{{1, 3}}", 0);
    assert_eval_eq("Length[First[FindIndependentVertexSet[CycleGraph[5]]]]", "2", 0);
    assert_eval_eq("FindIndependentVertexSet[StarGraph[5]]", "{{2, 3, 4, 5}}", 0);
    assert_eval_eq("Length[First[FindIndependentVertexSet[PathGraph[4]]]]", "2", 0);
    /* Duality: max independent set of g = max clique of the complement. */
    assert_eval_eq("Length[First[FindIndependentVertexSet[CycleGraph[5]]]] == Length[First[FindClique[GraphComplement[CycleGraph[5]]]]]", "True", 0);
    /* Direction ignored. */
    assert_eval_eq("Length[First[FindIndependentVertexSet[Graph[{1,2,3},{1->2,2->3}]]]]", "2", 0);
    assert_eval_eq("Head[FindIndependentVertexSet[5]]", "FindIndependentVertexSet", 0);
}

static void test_find_vertex_cover(void) {
    /* Edgeless → empty cover; complete K_n → n-1. */
    assert_eval_eq("FindVertexCover[Graph[{1,2,3},{}]]", "{}", 0);
    assert_eval_eq("Length[FindVertexCover[CompleteGraph[4]]]", "3", 0);
    /* Cycle / path / star / single edge sizes. */
    assert_eval_eq("FindVertexCover[CycleGraph[4]]", "{2, 4}", 0);
    assert_eval_eq("Length[FindVertexCover[PathGraph[4]]]", "2", 0);
    assert_eval_eq("FindVertexCover[StarGraph[5]]", "{1}", 0);
    assert_eval_eq("Length[FindVertexCover[Graph[{1,2},{1<->2}]]]", "1", 0);
    assert_eval_eq("Length[FindVertexCover[Graph[{1,2,3},{1<->2,2<->3,3<->1}]]]", "2", 0);
    /* Gallai identity: |min cover| + |max independent set| = n. */
    assert_eval_eq("Length[FindVertexCover[CycleGraph[5]]] + Length[First[FindIndependentVertexSet[CycleGraph[5]]]] == 5", "True", 0);
    assert_eval_eq("Length[FindVertexCover[CompleteGraph[5]]] + Length[First[FindIndependentVertexSet[CompleteGraph[5]]]] == 5", "True", 0);
    /* Direction ignored. */
    assert_eval_eq("Length[FindVertexCover[Graph[{1,2,3},{1->2,2->3}]]]", "1", 0);
    assert_eval_eq("Head[FindVertexCover[5]]", "FindVertexCover", 0);
}

static void test_graph_reciprocity(void) {
    /* Undirected graphs are fully reciprocal. */
    assert_eval_eq("GraphReciprocity[CycleGraph[4]]", "1", 0);
    assert_eval_eq("GraphReciprocity[CompleteGraph[5]]", "1", 0);
    assert_eval_eq("GraphReciprocity[Graph[{1,2},{1<->2}]]", "1", 0);
    /* Directed reciprocity is the fraction of reciprocated arcs. */
    assert_eval_eq("GraphReciprocity[Graph[{1,2,3},{1->2,2->3,3->1}]]", "0", 0);
    assert_eval_eq("GraphReciprocity[Graph[{1,2},{1->2,2->1}]]", "1", 0);
    assert_eval_eq("GraphReciprocity[Graph[{1,2,3},{1->2,2->1,2->3}]]", "2/3", 0);
    assert_eval_eq("GraphReciprocity[Graph[{1,2,3,4},{1->2,2->1,3->4}]]", "2/3", 0);
    assert_eval_eq("GraphReciprocity[Graph[{1,2},{1->2}]]", "0", 0);
    /* No edges → 0 by convention. */
    assert_eval_eq("GraphReciprocity[Graph[{1,2,3},{}]]", "0", 0);
    assert_eval_eq("Head[GraphReciprocity[5]]", "GraphReciprocity", 0);
}

static void test_chromatic_polynomial(void) {
    /* Numeric k gives the proper-colouring count. */
    assert_eval_eq("ChromaticPolynomial[Graph[{1,2,3},{}], 2]", "8", 0);   /* k^3 */
    assert_eval_eq("ChromaticPolynomial[CompleteGraph[3], 3]", "6", 0);    /* 3! */
    assert_eval_eq("ChromaticPolynomial[CompleteGraph[3], 2]", "0", 0);    /* not 2-colorable */
    assert_eval_eq("ChromaticPolynomial[CompleteGraph[4], 4]", "24", 0);
    assert_eval_eq("ChromaticPolynomial[PathGraph[3], 2]", "2", 0);
    assert_eval_eq("ChromaticPolynomial[CycleGraph[4], 3]", "18", 0);
    assert_eval_eq("ChromaticPolynomial[CycleGraph[5], 3]", "30", 0);
    assert_eval_eq("ChromaticPolynomial[Graph[{1},{}], 7]", "7", 0);
    /* Bipartite → 2-colorable; odd cycle → not. */
    assert_eval_eq("ChromaticPolynomial[CycleGraph[4], 2]", "2", 0);
    assert_eval_eq("ChromaticPolynomial[CycleGraph[5], 2]", "0", 0);
    /* Symbolic k yields a polynomial; substitute to check. */
    assert_eval_eq("ChromaticPolynomial[Graph[{1,2},{1<->2}], k] /. k->5", "20", 0);
    assert_eval_eq("ChromaticPolynomial[CompleteGraph[3], k] /. k->4", "24", 0);
    /* Oversized / non-graph stays unevaluated. */
    assert_eval_eq("Head[ChromaticPolynomial[CompleteGraph[9], k]]", "ChromaticPolynomial", 0);
    assert_eval_eq("Head[ChromaticPolynomial[5, k]]", "ChromaticPolynomial", 0);
}

static void test_chromatic_number(void) {
    assert_eval_eq("ChromaticNumber[Graph[{1,2,3},{}]]", "1", 0);
    assert_eval_eq("ChromaticNumber[Graph[{1},{}]]", "1", 0);
    assert_eval_eq("ChromaticNumber[CompleteGraph[4]]", "4", 0);
    assert_eval_eq("ChromaticNumber[CompleteGraph[5]]", "5", 0);
    assert_eval_eq("ChromaticNumber[CycleGraph[4]]", "2", 0);   /* bipartite */
    assert_eval_eq("ChromaticNumber[CycleGraph[5]]", "3", 0);   /* odd cycle */
    assert_eval_eq("ChromaticNumber[Graph[{1,2,3},{1<->2,2<->3,3<->1}]]", "3", 0);
    assert_eval_eq("ChromaticNumber[PathGraph[4]]", "2", 0);
    assert_eval_eq("ChromaticNumber[StarGraph[6]]", "2", 0);
    assert_eval_eq("ChromaticNumber[WheelGraph[6]]", "4", 0);   /* odd rim */
    assert_eval_eq("ChromaticNumber[WheelGraph[7]]", "3", 0);   /* even rim */
    /* Consistent with the chromatic polynomial: not 2-colorable ⇒ P(2)=0. */
    assert_eval_eq("ChromaticNumber[CycleGraph[5]] == 3 && ChromaticPolynomial[CycleGraph[5],2]==0", "True", 0);
    assert_eval_eq("ChromaticNumber[Graph[{1,2,3},{1->2,2->3,3->1}]]", "3", 0);
    assert_eval_eq("Head[ChromaticNumber[5]]", "ChromaticNumber", 0);
}

static void test_degree_sequence(void) {
    assert_eval_eq("DegreeSequence[CompleteGraph[4]]", "{3, 3, 3, 3}", 0);
    assert_eval_eq("DegreeSequence[StarGraph[5]]", "{4, 1, 1, 1, 1}", 0);
    assert_eval_eq("DegreeSequence[PathGraph[4]]", "{2, 2, 1, 1}", 0);
    assert_eval_eq("DegreeSequence[CycleGraph[5]]", "{2, 2, 2, 2, 2}", 0);
    assert_eval_eq("DegreeSequence[Graph[{1,2,3},{}]]", "{0, 0, 0}", 0);
    assert_eval_eq("DegreeSequence[Graph[{1,2,3,4},{1<->2,1<->3,1<->4,2<->3,2<->4}]]", "{3, 3, 2, 2}", 0);
    /* Directed: total degree, sorted descending. */
    assert_eval_eq("DegreeSequence[Graph[{1,2,3},{1->2,2->3}]]", "{2, 1, 1}", 0);
    /* Handshake lemma and permutation-of-DegreeCentrality consistency. */
    assert_eval_eq("Total[DegreeSequence[CompleteGraph[6]]] == 2*EdgeCount[CompleteGraph[6]]", "True", 0);
    assert_eval_eq("Sort[DegreeSequence[StarGraph[5]]] == Sort[DegreeCentrality[StarGraph[5]]]", "True", 0);
    assert_eval_eq("Head[DegreeSequence[5]]", "DegreeSequence", 0);
}

static void test_tree_graph_q(void) {
    assert_eval_eq("TreeGraphQ[PathGraph[4]]", "True", 0);
    assert_eval_eq("TreeGraphQ[StarGraph[5]]", "True", 0);
    assert_eval_eq("TreeGraphQ[Graph[{1},{}]]", "True", 0);
    assert_eval_eq("TreeGraphQ[Graph[{1,2,3,4,5},{1<->2,1<->3,2<->4,2<->5}]]", "True", 0);
    assert_eval_eq("TreeGraphQ[Graph[{1,2,3},{1->2,2->3}]]", "True", 0);
    /* Not trees: cycles, complete, edgeless multi-vertex, forests. */
    assert_eval_eq("TreeGraphQ[CycleGraph[4]]", "False", 0);
    assert_eval_eq("TreeGraphQ[CompleteGraph[4]]", "False", 0);
    assert_eval_eq("TreeGraphQ[Graph[{1,2,3},{}]]", "False", 0);
    assert_eval_eq("TreeGraphQ[Graph[{1,2,3,4},{1<->2,3<->4}]]", "False", 0);
    assert_eval_eq("TreeGraphQ[Graph[{1,2,3},{1<->2,2<->3,3<->1}]]", "False", 0);
    /* Characterisation: tree ⇔ connected and n-1 edges. */
    assert_eval_eq("TreeGraphQ[PathGraph[5]] == (ConnectedGraphQ[PathGraph[5]] && EdgeCount[PathGraph[5]]==4)", "True", 0);
    assert_eval_eq("Head[TreeGraphQ[5]]", "TreeGraphQ", 0);
}

static void test_strongly_connected_q(void) {
    /* Directed: needs cycles reaching both ways. */
    assert_eval_eq("StronglyConnectedGraphQ[Graph[{1,2,3},{1->2,2->3,3->1}]]", "True", 0);
    assert_eval_eq("StronglyConnectedGraphQ[Graph[{1,2,3},{1->2,2->3}]]", "False", 0);
    assert_eval_eq("StronglyConnectedGraphQ[Graph[{1,2},{1->2,2->1}]]", "True", 0);
    assert_eval_eq("StronglyConnectedGraphQ[Graph[{1,2,3},{1->2,2->1,2->3}]]", "False", 0);
    /* Undirected: coincides with connectivity. */
    assert_eval_eq("StronglyConnectedGraphQ[CycleGraph[4]]", "True", 0);
    assert_eval_eq("StronglyConnectedGraphQ[PathGraph[4]]", "True", 0);
    assert_eval_eq("StronglyConnectedGraphQ[Graph[{1,2,3},{1<->2}]]", "False", 0);
    assert_eval_eq("StronglyConnectedGraphQ[Graph[{1},{}]]", "True", 0);
    assert_eval_eq("StronglyConnectedGraphQ[Graph[{1,2},{}]]", "False", 0);
    /* Consistent with a single all-covering strongly connected component. */
    assert_eval_eq("StronglyConnectedGraphQ[Graph[{1,2,3},{1->2,2->3,3->1}]] == (Length[StronglyConnectedComponents[Graph[{1,2,3},{1->2,2->3,3->1}]]] == 1)", "True", 0);
    assert_eval_eq("Head[StronglyConnectedGraphQ[5]]", "StronglyConnectedGraphQ", 0);
}

static void test_hamiltonian_graph_q(void) {
    assert_eval_eq("HamiltonianGraphQ[CycleGraph[5]]", "True", 0);
    assert_eval_eq("HamiltonianGraphQ[CompleteGraph[4]]", "True", 0);
    assert_eval_eq("HamiltonianGraphQ[WheelGraph[5]]", "True", 0);
    assert_eval_eq("HamiltonianGraphQ[Graph[{1,2,3},{1->2,2->3,3->1}]]", "True", 0);
    /* Not Hamiltonian: paths, stars, broken directed, edgeless, n<3. */
    assert_eval_eq("HamiltonianGraphQ[PathGraph[4]]", "False", 0);
    assert_eval_eq("HamiltonianGraphQ[StarGraph[5]]", "False", 0);
    assert_eval_eq("HamiltonianGraphQ[Graph[{1,2,3},{1->2,2->3,1->3}]]", "False", 0);
    assert_eval_eq("HamiltonianGraphQ[Graph[{1,2,3},{}]]", "False", 0);
    assert_eval_eq("HamiltonianGraphQ[Graph[{1,2},{1<->2}]]", "False", 0);
    /* Agrees with FindHamiltonianCycle. */
    assert_eval_eq("HamiltonianGraphQ[CycleGraph[5]] == (FindHamiltonianCycle[CycleGraph[5]] =!= {})", "True", 0);
    assert_eval_eq("HamiltonianGraphQ[PathGraph[4]] == (FindHamiltonianCycle[PathGraph[4]] =!= {})", "True", 0);
    assert_eval_eq("Head[HamiltonianGraphQ[5]]", "HamiltonianGraphQ", 0);
}

static void test_regular_graph_q(void) {
    assert_eval_eq("RegularGraphQ[CycleGraph[5]]", "True", 0);
    assert_eval_eq("RegularGraphQ[CompleteGraph[4]]", "True", 0);
    assert_eval_eq("RegularGraphQ[Graph[{1,2,3},{}]]", "True", 0);   /* 0-regular */
    assert_eval_eq("RegularGraphQ[Graph[{1},{}]]", "True", 0);
    assert_eval_eq("RegularGraphQ[Graph[{1,2,3,4,5,6},{1<->4,1<->5,1<->6,2<->4,2<->5,2<->6,3<->4,3<->5,3<->6}]]", "True", 0);
    /* Not regular. */
    assert_eval_eq("RegularGraphQ[PathGraph[4]]", "False", 0);
    assert_eval_eq("RegularGraphQ[StarGraph[5]]", "False", 0);
    assert_eval_eq("RegularGraphQ[WheelGraph[5]]", "False", 0);
    /* Directed: equal in- and out-degrees. */
    assert_eval_eq("RegularGraphQ[Graph[{1,2,3},{1->2,2->3,3->1}]]", "True", 0);
    assert_eval_eq("RegularGraphQ[Graph[{1,2,3},{1->2,2->3}]]", "False", 0);
    assert_eval_eq("Head[RegularGraphQ[5]]", "RegularGraphQ", 0);
}

static void test_complete_graph_q(void) {
    assert_eval_eq("CompleteGraphQ[CompleteGraph[4]]", "True", 0);
    assert_eval_eq("CompleteGraphQ[CompleteGraph[5]]", "True", 0);
    assert_eval_eq("CompleteGraphQ[CycleGraph[3]]", "True", 0);   /* C3 = K3 */
    assert_eval_eq("CompleteGraphQ[Graph[{1,2},{1<->2}]]", "True", 0);
    assert_eval_eq("CompleteGraphQ[Graph[{1},{}]]", "True", 0);
    assert_eval_eq("CompleteGraphQ[Graph[{1,2,3},{1->2,2->3,3->1}]]", "True", 0);
    /* Not complete. */
    assert_eval_eq("CompleteGraphQ[CycleGraph[4]]", "False", 0);
    assert_eval_eq("CompleteGraphQ[PathGraph[4]]", "False", 0);
    assert_eval_eq("CompleteGraphQ[Graph[{1,2,3},{}]]", "False", 0);
    assert_eval_eq("CompleteGraphQ[Graph[{1,2,3,4},{1<->2,1<->3,1<->4,2<->3,2<->4}]]", "False", 0);
    assert_eval_eq("Head[CompleteGraphQ[5]]", "CompleteGraphQ", 0);
}

static void test_graph_union(void) {
    /* Union merges vertices and edges. */
    assert_eval_eq("VertexList[GraphUnion[Graph[{1,2},{1<->2}], Graph[{2,3},{2<->3}]]]", "{1, 2, 3}", 0);
    assert_eval_eq("VertexCount[GraphUnion[Graph[{1,2,3},{1<->2,2<->3}], Graph[{3,4,5},{3<->4,4<->5}]]]", "5", 0);
    assert_eval_eq("EdgeCount[GraphUnion[Graph[{1,2,3},{1<->2,2<->3}], Graph[{3,4,5},{3<->4,4<->5}]]]", "4", 0);
    /* Shared / symmetric edges deduped; self-union is idempotent. */
    assert_eval_eq("EdgeCount[GraphUnion[Graph[{1,2,3},{1<->2,2<->3}], Graph[{2,3,4},{2<->3,3<->4}]]]", "3", 0);
    assert_eval_eq("EdgeCount[GraphUnion[Graph[{1,2},{1<->2}], Graph[{1,2},{2<->1}]]]", "1", 0);
    assert_eval_eq("EdgeCount[GraphUnion[CycleGraph[4], CycleGraph[4]]] == EdgeCount[CycleGraph[4]]", "True", 0);
    /* Directed edges are not symmetric. */
    assert_eval_eq("EdgeCount[GraphUnion[Graph[{1,2},{1->2}], Graph[{1,2},{2->1}]]]", "2", 0);
    /* Edgeless union and validity. */
    assert_eval_eq("VertexCount[GraphUnion[Graph[{1,2},{}], Graph[{3,4},{}]]]", "4", 0);
    assert_eval_eq("GraphQ[GraphUnion[PathGraph[3], CycleGraph[3]]]", "True", 0);
    assert_eval_eq("Head[GraphUnion[5, CycleGraph[3]]]", "GraphUnion", 0);
}

static void test_graph_intersection(void) {
    /* Common vertices and common edges. */
    assert_eval_eq("VertexList[GraphIntersection[Graph[{1,2,3},{1<->2,2<->3}], Graph[{2,3,4},{2<->3,3<->4}]]]", "{2, 3}", 0);
    assert_eval_eq("EdgeList[GraphIntersection[Graph[{1,2,3},{1<->2,2<->3}], Graph[{2,3,4},{2<->3,3<->4}]]]", "{2 <-> 3}", 0);
    /* Identical graphs → self. */
    assert_eval_eq("EdgeCount[GraphIntersection[CycleGraph[4], CycleGraph[4]]]", "4", 0);
    assert_eval_eq("VertexCount[GraphIntersection[CycleGraph[4], CycleGraph[4]]]", "4", 0);
    /* Disjoint / no shared edges. */
    assert_eval_eq("VertexCount[GraphIntersection[Graph[{1,2},{1<->2}], Graph[{3,4},{3<->4}]]]", "0", 0);
    assert_eval_eq("EdgeCount[GraphIntersection[Graph[{1,2,3},{1<->2}], Graph[{1,2,3},{2<->3}]]]", "0", 0);
    /* Symmetric undirected common; directed not symmetric. */
    assert_eval_eq("EdgeCount[GraphIntersection[Graph[{1,2},{1<->2}], Graph[{1,2},{2<->1}]]]", "1", 0);
    assert_eval_eq("EdgeCount[GraphIntersection[Graph[{1,2},{1->2}], Graph[{1,2},{2->1}]]]", "0", 0);
    assert_eval_eq("EdgeCount[GraphIntersection[CompleteGraph[4], CycleGraph[4]]]", "4", 0);
    assert_eval_eq("Head[GraphIntersection[5, CycleGraph[3]]]", "GraphIntersection", 0);
}

static void test_graph_difference(void) {
    /* K4 minus C4 leaves the two diagonals, keeps all vertices. */
    assert_eval_eq("EdgeCount[GraphDifference[CompleteGraph[4], CycleGraph[4]]]", "2", 0);
    assert_eval_eq("VertexCount[GraphDifference[CompleteGraph[4], CycleGraph[4]]]", "4", 0);
    /* Self-difference: edgeless on the same vertices. */
    assert_eval_eq("EdgeCount[GraphDifference[CycleGraph[5], CycleGraph[5]]]", "0", 0);
    assert_eval_eq("VertexCount[GraphDifference[CycleGraph[5], CycleGraph[5]]]", "5", 0);
    /* Disjoint g2 removes nothing; a shared edge is removed. */
    assert_eval_eq("EdgeCount[GraphDifference[Graph[{1,2,3},{1<->2,2<->3}], Graph[{4,5},{4<->5}]]]", "2", 0);
    assert_eval_eq("EdgeList[GraphDifference[Graph[{1,2,3},{1<->2,2<->3}], Graph[{2,3},{2<->3}]]]", "{1 <-> 2}", 0);
    /* Symmetric undirected removal; directed not symmetric. */
    assert_eval_eq("EdgeCount[GraphDifference[Graph[{1,2},{1<->2}], Graph[{1,2},{2<->1}]]]", "0", 0);
    assert_eval_eq("EdgeCount[GraphDifference[Graph[{1,2},{1->2}], Graph[{1,2},{2->1}]]]", "1", 0);
    assert_eval_eq("GraphQ[GraphDifference[CompleteGraph[5], CycleGraph[5]]]", "True", 0);
    assert_eval_eq("Head[GraphDifference[5, CycleGraph[3]]]", "GraphDifference", 0);
}

static void test_graph_reverse(void) {
    assert_eval_eq("EdgeList[ReverseGraph[Graph[{1,2,3},{1->2,2->3}]]]", "{2 -> 1, 3 -> 2}", 0);
    assert_eval_eq("VertexList[ReverseGraph[Graph[{1,2,3},{1->2,2->3}]]]", "{1, 2, 3}", 0);
    assert_eval_eq("EdgeCount[ReverseGraph[Graph[{1,2,3},{1->2,2->3,3->1}]]]", "3", 0);
    /* Undirected edges are unchanged; reversal is an involution on directed graphs. */
    assert_eval_eq("EdgeList[ReverseGraph[Graph[{1,2},{1<->2}]]]", "{1 <-> 2}", 0);
    assert_eval_eq("ReverseGraph[ReverseGraph[Graph[{1,2,3},{1->2,2->3}]]] === Graph[{1,2,3},{1->2,2->3}]", "True", 0);
    /* Swaps in- and out-degree. */
    assert_eval_eq("VertexOutDegree[ReverseGraph[Graph[{1,2,3},{1->2,1->3}]], 1]", "0", 0);
    assert_eval_eq("VertexInDegree[ReverseGraph[Graph[{1,2,3},{1->2,1->3}]], 1]", "2", 0);
    /* A reversed directed cycle is still strongly connected. */
    assert_eval_eq("StronglyConnectedGraphQ[ReverseGraph[Graph[{1,2,3},{1->2,2->3,3->1}]]]", "True", 0);
    assert_eval_eq("GraphQ[ReverseGraph[Graph[{1,2},{1->2}]]]", "True", 0);
    assert_eval_eq("Head[ReverseGraph[5]]", "ReverseGraph", 0);
}

static void test_path_graph_q(void) {
    assert_eval_eq("PathGraphQ[PathGraph[4]]", "True", 0);
    assert_eval_eq("PathGraphQ[Graph[{1},{}]]", "True", 0);
    assert_eval_eq("PathGraphQ[Graph[{1,2},{1<->2}]]", "True", 0);
    assert_eval_eq("PathGraphQ[Graph[{1,2,3},{1->2,2->3}]]", "True", 0);
    /* Not paths: cycles, stars, branches, disconnected, edgeless. */
    assert_eval_eq("PathGraphQ[CycleGraph[4]]", "False", 0);
    assert_eval_eq("PathGraphQ[StarGraph[4]]", "False", 0);
    assert_eval_eq("PathGraphQ[CompleteGraph[3]]", "False", 0);
    assert_eval_eq("PathGraphQ[Graph[{1,2,3,4},{1<->2,2<->3,2<->4}]]", "False", 0);
    assert_eval_eq("PathGraphQ[Graph[{1,2,3,4},{1<->2,3<->4}]]", "False", 0);
    assert_eval_eq("PathGraphQ[Graph[{1,2,3},{}]]", "False", 0);
    /* A path is a tree. */
    assert_eval_eq("PathGraphQ[PathGraph[5]] && TreeGraphQ[PathGraph[5]]", "True", 0);
    assert_eval_eq("Head[PathGraphQ[5]]", "PathGraphQ", 0);
}

static void test_vertex_contract(void) {
    /* Contracting an edge of a triangle yields a single edge. */
    assert_eval_eq("EdgeList[VertexContract[Graph[{1,2,3},{1<->2,2<->3,3<->1}], {1,2}]]", "{1 <-> 3}", 0);
    assert_eval_eq("VertexList[VertexContract[Graph[{1,2,3},{1<->2,2<->3,3<->1}], {1,2}]]", "{1, 3}", 0);
    /* Singleton contraction is a no-op. */
    assert_eval_eq("EdgeCount[VertexContract[CycleGraph[4], {1}]]", "4", 0);
    assert_eval_eq("VertexCount[VertexContract[CycleGraph[4], {1}]]", "4", 0);
    /* Contracting a path's endpoints closes it into a triangle. */
    assert_eval_eq("EdgeCount[VertexContract[PathGraph[4], {1,4}]]", "3", 0);
    assert_eval_eq("VertexCount[VertexContract[PathGraph[4], {1,4}]]", "3", 0);
    /* Contract all → one vertex, no edges. */
    assert_eval_eq("VertexCount[VertexContract[CompleteGraph[4], {1,2,3,4}]]", "1", 0);
    assert_eval_eq("EdgeCount[VertexContract[CompleteGraph[4], {1,2,3,4}]]", "0", 0);
    /* Parallel edges collapse; directed self-loops drop. */
    assert_eval_eq("EdgeCount[VertexContract[Graph[{1,2,3},{1<->3,2<->3}], {1,2}]]", "1", 0);
    assert_eval_eq("EdgeCount[VertexContract[Graph[{1,2,3},{1->2,2->3}], {1,2}]]", "1", 0);
    assert_eval_eq("GraphQ[VertexContract[CompleteGraph[4], {1,2}]]", "True", 0);
    assert_eval_eq("Head[VertexContract[CycleGraph[3], {9}]]", "VertexContract", 0);
    assert_eval_eq("Head[VertexContract[5, {1}]]", "VertexContract", 0);
}

static void test_pagerank_centrality(void) {
    /* Exact rational vector summing to 1. */
    assert_eval_eq("Total[PageRankCentrality[CycleGraph[4]]]", "1", 0);
    assert_eval_eq("Total[PageRankCentrality[StarGraph[5]]]", "1", 0);
    assert_eval_eq("Total[PageRankCentrality[Graph[{1,2,3},{1->2,2->3,3->1}]]]", "1", 0);
    /* Regular graphs / all-dangling → uniform 1/n. */
    assert_eval_eq("PageRankCentrality[CycleGraph[4]]", "{1/4, 1/4, 1/4, 1/4}", 0);
    assert_eval_eq("PageRankCentrality[CompleteGraph[5]]", "{1/5, 1/5, 1/5, 1/5, 1/5}", 0);
    assert_eval_eq("PageRankCentrality[Graph[{1,2,3},{}]]", "{1/3, 1/3, 1/3}", 0);
    assert_eval_eq("PageRankCentrality[Graph[{1},{}]]", "{1}", 0);
    /* Star: exact rationals, centre outranks the leaves. */
    assert_eval_eq("PageRankCentrality[StarGraph[4]]", "{71/148, 77/444, 77/444, 77/444}", 0);
    assert_eval_eq("First[PageRankCentrality[StarGraph[5]]] > PageRankCentrality[StarGraph[5]][[2]]", "True", 0);
    /* Directed hub. */
    assert_eval_eq("PageRankCentrality[Graph[{1,2,3},{1->2,1->3}]]", "{20/77, 57/154, 57/154}", 0);
    assert_eval_eq("Head[PageRankCentrality[5]]", "PageRankCentrality", 0);
}

static void test_katz_centrality(void) {
    /* alpha = 0 → base weights only. */
    assert_eval_eq("KatzCentrality[CycleGraph[4], 0]", "{1, 1, 1, 1}", 0);
    assert_eval_eq("KatzCentrality[Graph[{1,2,3},{}], 1/5]", "{1, 1, 1}", 0);
    assert_eval_eq("KatzCentrality[Graph[{1},{}], 1/2]", "{1}", 0);
    /* Exact rationals; regular graph is uniform, path centre outranks ends. */
    assert_eval_eq("KatzCentrality[Graph[{1,2},{1<->2}], 1/10]", "{10/9, 10/9}", 0);
    assert_eval_eq("KatzCentrality[CycleGraph[4], 1/10]", "{5/4, 5/4, 5/4, 5/4}", 0);
    assert_eval_eq("KatzCentrality[PathGraph[3], 1/10]", "{55/49, 60/49, 55/49}", 0);
    /* Directed: score comes from in-edges (who points at you). */
    assert_eval_eq("KatzCentrality[Graph[{1,2,3},{1->2,1->3}], 1/10][[2]] > KatzCentrality[Graph[{1,2,3},{1->2,1->3}], 1/10][[1]]", "True", 0);
    /* Non-numeric alpha / non-graph / arity stay unevaluated. */
    assert_eval_eq("Head[KatzCentrality[CycleGraph[3], a]]", "KatzCentrality", 0);
    assert_eval_eq("Head[KatzCentrality[5, 1/10]]", "KatzCentrality", 0);
    assert_eval_eq("Head[KatzCentrality[CycleGraph[3]]]", "KatzCentrality", 0);
}

static void test_graph_join(void) {
    /* K1 join K1 = a single edge. */
    assert_eval_eq("EdgeCount[GraphJoin[Graph[{1},{}], Graph[{1},{}]]]", "1", 0);
    assert_eval_eq("VertexCount[GraphJoin[Graph[{1},{}], Graph[{1},{}]]]", "2", 0);
    /* Edge count = m1 + m2 + n1*n2; vertices relabelled 1..n1+n2. */
    assert_eval_eq("EdgeCount[GraphJoin[PathGraph[3], PathGraph[2]]]", "9", 0);
    assert_eval_eq("VertexCount[GraphJoin[PathGraph[3], PathGraph[2]]]", "5", 0);
    assert_eval_eq("VertexList[GraphJoin[Graph[{a,b},{a<->b}], Graph[{x},{}]]]", "{1, 2, 3}", 0);
    /* P2 join K1 is a triangle; Km join Kn is K(m+n). */
    assert_eval_eq("CompleteGraphQ[GraphJoin[PathGraph[2], Graph[{1},{}]]]", "True", 0);
    assert_eval_eq("CompleteGraphQ[GraphJoin[CompleteGraph[2], CompleteGraph[2]]]", "True", 0);
    assert_eval_eq("CompleteGraphQ[GraphJoin[CompleteGraph[2], CompleteGraph[3]]]", "True", 0);
    assert_eval_eq("GraphQ[GraphJoin[CycleGraph[3], CycleGraph[3]]]", "True", 0);
    assert_eval_eq("Head[GraphJoin[5, CycleGraph[3]]]", "GraphJoin", 0);
}

static void test_index_graph(void) {
    /* Symbolic labels → 1..n, edges remapped, kinds preserved. */
    assert_eval_eq("VertexList[IndexGraph[Graph[{x,y,z},{x<->z}]]]", "{1, 2, 3}", 0);
    assert_eval_eq("EdgeList[IndexGraph[Graph[{a,b,c},{a<->b,b<->c}]]]", "{1 <-> 2, 2 <-> 3}", 0);
    assert_eval_eq("EdgeList[IndexGraph[Graph[{a,b,c},{a->b,b->c}]]]", "{1 -> 2, 2 -> 3}", 0);
    assert_eval_eq("VertexList[IndexGraph[Graph[{a,b},{a<->b}], 0]]", "{0, 1}", 0);
    assert_eval_eq("EdgeCount[IndexGraph[CompleteGraph[4]]]", "6", 0);
    /* Already 1..n is unchanged; structure (edge count) is preserved. */
    assert_eval_eq("IndexGraph[PathGraph[3]] === PathGraph[3]", "True", 0);
    assert_eval_eq("EdgeCount[IndexGraph[Graph[{p,q,r},{p<->q,q<->r,r<->p}]]]", "3", 0);
    assert_eval_eq("GraphQ[IndexGraph[Graph[{a,b,c},{a<->b}]]]", "True", 0);
    /* Non-integer start / non-graph stay unevaluated. */
    assert_eval_eq("Head[IndexGraph[CycleGraph[3], x]]", "IndexGraph", 0);
    assert_eval_eq("Head[IndexGraph[5]]", "IndexGraph", 0);
}

static void test_empty_and_mixed_q(void) {
    /* EmptyGraphQ: no edges. */
    assert_eval_eq("EmptyGraphQ[Graph[{1,2,3},{}]]", "True", 0);
    assert_eval_eq("EmptyGraphQ[Graph[{1},{}]]", "True", 0);
    assert_eval_eq("EmptyGraphQ[CycleGraph[3]]", "False", 0);
    assert_eval_eq("EmptyGraphQ[Graph[{1,2},{1<->2}]]", "False", 0);
    assert_eval_eq("Head[EmptyGraphQ[5]]", "EmptyGraphQ", 0);
    /* MixedGraphQ: both directed and undirected edges present. */
    assert_eval_eq("MixedGraphQ[Graph[{1,2,3},{1->2,2<->3}]]", "True", 0);
    assert_eval_eq("MixedGraphQ[Graph[{1,2,3},{1->2,2->3}]]", "False", 0);
    assert_eval_eq("MixedGraphQ[CycleGraph[4]]", "False", 0);
    assert_eval_eq("MixedGraphQ[Graph[{1,2},{}]]", "False", 0);
    assert_eval_eq("Head[MixedGraphQ[5]]", "MixedGraphQ", 0);
}

static void test_graph_product(void) {
    /* P2 [] P2 = C4. */
    assert_eval_eq("VertexCount[GraphProduct[Graph[{1,2},{1<->2}], Graph[{1,2},{1<->2}], \"Cartesian\"]]", "4", 0);
    assert_eval_eq("EdgeCount[GraphProduct[Graph[{1,2},{1<->2}], Graph[{1,2},{1<->2}], \"Cartesian\"]]", "4", 0);
    /* Tensor / Strong / Lexicographic of K2 with K2. */
    assert_eval_eq("EdgeCount[GraphProduct[Graph[{1,2},{1<->2}], Graph[{1,2},{1<->2}], \"Tensor\"]]", "2", 0);
    assert_eval_eq("CompleteGraphQ[GraphProduct[Graph[{1,2},{1<->2}], Graph[{1,2},{1<->2}], \"Strong\"]]", "True", 0);
    assert_eval_eq("EdgeCount[GraphProduct[Graph[{1,2},{1<->2}], Graph[{1,2},{1<->2}], \"Lexicographic\"]]", "6", 0);
    /* C4 [] K2 is the 3-regular cube (8 vertices, 12 edges). */
    assert_eval_eq("VertexCount[GraphProduct[CycleGraph[4], Graph[{1,2},{1<->2}], \"Cartesian\"]]", "8", 0);
    assert_eval_eq("EdgeCount[GraphProduct[CycleGraph[4], Graph[{1,2},{1<->2}], \"Cartesian\"]]", "12", 0);
    assert_eval_eq("RegularGraphQ[GraphProduct[CycleGraph[4], Graph[{1,2},{1<->2}], \"Cartesian\"]]", "True", 0);
    /* Cartesian grid P3 x P2 has 3*1 + 2*2 = 7 edges. */
    assert_eval_eq("EdgeCount[GraphProduct[PathGraph[3], Graph[{1,2},{1<->2}], \"Cartesian\"]]", "7", 0);
    /* Unknown type / non-string / non-graph stay unevaluated. */
    assert_eval_eq("Head[GraphProduct[Graph[{1,2},{1<->2}], Graph[{1,2},{1<->2}], \"Nope\"]]", "GraphProduct", 0);
    assert_eval_eq("Head[GraphProduct[5, Graph[{1,2},{1<->2}], \"Cartesian\"]]", "GraphProduct", 0);
}

static void test_turan_graph(void) {
    assert_eval_eq("EdgeCount[TuranGraph[4,2]]", "4", 0);    /* C4 */
    assert_eval_eq("EdgeCount[TuranGraph[5,2]]", "6", 0);    /* K_{2,3} */
    assert_eval_eq("EdgeCount[TuranGraph[6,3]]", "12", 0);   /* octahedron */
    assert_eval_eq("VertexCount[TuranGraph[7,3]]", "7", 0);
    assert_eval_eq("EdgeCount[TuranGraph[5,1]]", "0", 0);    /* edgeless */
    assert_eval_eq("CompleteGraphQ[TuranGraph[5,5]]", "True", 0);   /* K_n */
    assert_eval_eq("RegularGraphQ[TuranGraph[6,3]]", "True", 0);
    /* Cross-checks with other builtins. */
    assert_eval_eq("BipartiteGraphQ[TuranGraph[4,2]]", "True", 0);
    assert_eval_eq("ChromaticNumber[TuranGraph[4,2]]", "2", 0);
    assert_eval_eq("ChromaticNumber[TuranGraph[6,3]]", "3", 0);
    /* Bad arguments stay unevaluated. */
    assert_eval_eq("Head[TuranGraph[5,0]]", "TuranGraph", 0);
    assert_eval_eq("Head[TuranGraph[5,x]]", "TuranGraph", 0);
}

static void test_complete_kary_tree(void) {
    /* Binary trees (default k=2). */
    assert_eval_eq("VertexCount[CompleteKaryTree[2]]", "3", 0);
    assert_eval_eq("VertexCount[CompleteKaryTree[3]]", "7", 0);
    assert_eval_eq("EdgeCount[CompleteKaryTree[3]]", "6", 0);
    assert_eval_eq("EdgeList[CompleteKaryTree[2]]", "{1 <-> 2, 1 <-> 3}", 0);
    assert_eval_eq("TreeGraphQ[CompleteKaryTree[3]]", "True", 0);
    /* k-ary: L3 ternary = 13 vertices; root has k children. */
    assert_eval_eq("VertexCount[CompleteKaryTree[3,3]]", "13", 0);
    assert_eval_eq("VertexOutDegree[CompleteKaryTree[2,3], 1]", "3", 0);
    /* k=1 is a path; L=1 is a single vertex. */
    assert_eval_eq("PathGraphQ[CompleteKaryTree[5,1]]", "True", 0);
    assert_eval_eq("VertexCount[CompleteKaryTree[1]]", "1", 0);
    assert_eval_eq("EdgeCount[CompleteKaryTree[1]]", "0", 0);
    /* Bad arguments stay unevaluated. */
    assert_eval_eq("Head[CompleteKaryTree[0]]", "CompleteKaryTree", 0);
    assert_eval_eq("Head[CompleteKaryTree[x]]", "CompleteKaryTree", 0);
}

static void test_circulant_graph(void) {
    /* C_n({1}) is the cycle; single-offset integer form too. */
    assert_eval_eq("EdgeCount[CirculantGraph[5,{1}]]", "5", 0);
    assert_eval_eq("EdgeCount[CirculantGraph[6,1]]", "6", 0);
    assert_eval_eq("EdgeCount[CirculantGraph[7,{1}]] == EdgeCount[CycleGraph[7]]", "True", 0);
    assert_eval_eq("RegularGraphQ[CirculantGraph[5,{1}]]", "True", 0);
    /* Full jump set → complete; other regular cases. */
    assert_eval_eq("CompleteGraphQ[CirculantGraph[6,{1,2,3}]]", "True", 0);
    assert_eval_eq("RegularGraphQ[CirculantGraph[8,{1,2}]]", "True", 0);
    assert_eval_eq("EdgeCount[CirculantGraph[8,{1,2}]]", "16", 0);
    /* A jump of n/2 contributes a single matching edge per vertex. */
    assert_eval_eq("EdgeCount[CirculantGraph[6,{3}]]", "3", 0);
    assert_eval_eq("VertexCount[CirculantGraph[7,{1,2}]]", "7", 0);
    /* Bad arguments stay unevaluated. */
    assert_eval_eq("Head[CirculantGraph[0,{1}]]", "CirculantGraph", 0);
    assert_eval_eq("Head[CirculantGraph[5,{x}]]", "CirculantGraph", 0);
}

static void test_ladder_graph(void) {
    assert_eval_eq("EdgeCount[LadderGraph[1]]", "1", 0);      /* K2 */
    assert_eval_eq("VertexCount[LadderGraph[1]]", "2", 0);
    assert_eval_eq("EdgeCount[LadderGraph[2]]", "4", 0);      /* C4 */
    assert_eval_eq("VertexCount[LadderGraph[3]]", "6", 0);
    assert_eval_eq("EdgeCount[LadderGraph[3]]", "7", 0);
    assert_eval_eq("EdgeCount[LadderGraph[5]] == 3*5-2", "True", 0);
    assert_eval_eq("ConnectedGraphQ[LadderGraph[4]]", "True", 0);
    assert_eval_eq("BipartiteGraphQ[LadderGraph[4]]", "True", 0);
    assert_eval_eq("TreeGraphQ[LadderGraph[3]]", "False", 0);
    /* Matches the Cartesian product P_n x P_2. */
    assert_eval_eq("EdgeCount[LadderGraph[3]] == EdgeCount[GraphProduct[PathGraph[3], Graph[{1,2},{1<->2}], \"Cartesian\"]]", "True", 0);
    assert_eval_eq("Head[LadderGraph[0]]", "LadderGraph", 0);
    assert_eval_eq("Head[LadderGraph[x]]", "LadderGraph", 0);
}

static void test_cocktail_party_graph(void) {
    assert_eval_eq("EdgeCount[CocktailPartyGraph[2]]", "4", 0);     /* C4 */
    assert_eval_eq("EdgeCount[CocktailPartyGraph[3]]", "12", 0);    /* octahedron */
    assert_eval_eq("VertexCount[CocktailPartyGraph[3]]", "6", 0);
    assert_eval_eq("EdgeCount[CocktailPartyGraph[1]]", "0", 0);
    assert_eval_eq("EdgeCount[CocktailPartyGraph[4]] == 2*4*3", "True", 0);
    assert_eval_eq("RegularGraphQ[CocktailPartyGraph[4]]", "True", 0);
    assert_eval_eq("First[DegreeCentrality[CocktailPartyGraph[4]]]", "6", 0);  /* 2n-2 */
    assert_eval_eq("ChromaticNumber[CocktailPartyGraph[3]]", "3", 0);
    /* Same as the balanced complete n-partite Turán graph on 2n vertices. */
    assert_eval_eq("EdgeCount[CocktailPartyGraph[3]] == EdgeCount[TuranGraph[6,3]]", "True", 0);
    assert_eval_eq("Head[CocktailPartyGraph[0]]", "CocktailPartyGraph", 0);
    assert_eval_eq("Head[CocktailPartyGraph[x]]", "CocktailPartyGraph", 0);
}

static void test_kneser_graph(void) {
    /* K(5,2) is the Petersen graph. */
    assert_eval_eq("VertexCount[KneserGraph[5,2]]", "10", 0);
    assert_eval_eq("EdgeCount[KneserGraph[5,2]]", "15", 0);
    assert_eval_eq("RegularGraphQ[KneserGraph[5,2]]", "True", 0);
    assert_eval_eq("ConnectedGraphQ[KneserGraph[5,2]]", "True", 0);
    assert_eval_eq("BipartiteGraphQ[KneserGraph[5,2]]", "False", 0);
    /* K(n,1) = K_n; K(4,2) = perfect matching; vertex labels are subsets. */
    assert_eval_eq("CompleteGraphQ[KneserGraph[5,1]]", "True", 0);
    assert_eval_eq("EdgeCount[KneserGraph[4,2]]", "3", 0);
    assert_eval_eq("RegularGraphQ[KneserGraph[4,2]]", "True", 0);
    assert_eval_eq("First[VertexList[KneserGraph[4,2]]]", "{1, 2}", 0);
    assert_eval_eq("VertexCount[KneserGraph[3,0]]", "1", 0);
    /* Bad arguments stay unevaluated. */
    assert_eval_eq("Head[KneserGraph[2,3]]", "KneserGraph", 0);
    assert_eval_eq("Head[KneserGraph[5,x]]", "KneserGraph", 0);
}

static void test_generalized_petersen_graph(void) {
    /* GP(5,2) is the Petersen graph. */
    assert_eval_eq("VertexCount[GeneralizedPetersenGraph[5,2]]", "10", 0);
    assert_eval_eq("EdgeCount[GeneralizedPetersenGraph[5,2]]", "15", 0);
    assert_eval_eq("RegularGraphQ[GeneralizedPetersenGraph[5,2]]", "True", 0);
    assert_eval_eq("BipartiteGraphQ[GeneralizedPetersenGraph[5,2]]", "False", 0);
    assert_eval_eq("EdgeCount[GeneralizedPetersenGraph[5,2]] == EdgeCount[KneserGraph[5,2]]", "True", 0);
    /* GP(4,1) is the cube; GP(n,1) is the n-prism (3n edges). */
    assert_eval_eq("VertexCount[GeneralizedPetersenGraph[4,1]]", "8", 0);
    assert_eval_eq("EdgeCount[GeneralizedPetersenGraph[4,1]]", "12", 0);
    assert_eval_eq("BipartiteGraphQ[GeneralizedPetersenGraph[4,1]]", "True", 0);
    assert_eval_eq("EdgeCount[GeneralizedPetersenGraph[6,1]]", "18", 0);
    assert_eval_eq("RegularGraphQ[GeneralizedPetersenGraph[3,1]]", "True", 0);
    assert_eval_eq("ConnectedGraphQ[GeneralizedPetersenGraph[8,3]]", "True", 0);
    /* Bad arguments stay unevaluated. */
    assert_eval_eq("Head[GeneralizedPetersenGraph[2,1]]", "GeneralizedPetersenGraph", 0);
    assert_eval_eq("Head[GeneralizedPetersenGraph[5,5]]", "GeneralizedPetersenGraph", 0);
}

static void test_friendship_graph(void) {
    assert_eval_eq("CompleteGraphQ[FriendshipGraph[1]]", "True", 0);   /* K3 */
    assert_eval_eq("VertexCount[FriendshipGraph[1]]", "3", 0);
    assert_eval_eq("VertexCount[FriendshipGraph[2]]", "5", 0);         /* bowtie */
    assert_eval_eq("EdgeCount[FriendshipGraph[2]]", "6", 0);
    assert_eval_eq("VertexCount[FriendshipGraph[4]]", "9", 0);         /* 2n+1 */
    assert_eval_eq("EdgeCount[FriendshipGraph[4]]", "12", 0);          /* 3n */
    assert_eval_eq("First[DegreeCentrality[FriendshipGraph[3]]]", "6", 0);  /* hub = 2n */
    assert_eval_eq("DegreeCentrality[FriendshipGraph[3]][[2]]", "2", 0);
    assert_eval_eq("ConnectedGraphQ[FriendshipGraph[3]]", "True", 0);
    assert_eval_eq("BipartiteGraphQ[FriendshipGraph[2]]", "False", 0);
    assert_eval_eq("ChromaticNumber[FriendshipGraph[3]]", "3", 0);
    assert_eval_eq("Head[FriendshipGraph[0]]", "FriendshipGraph", 0);
}

static void test_vertex_coreness(void) {
    assert_eval_eq("VertexCoreness[CompleteGraph[4]]", "{3, 3, 3, 3}", 0);
    assert_eval_eq("VertexCoreness[CycleGraph[4]]", "{2, 2, 2, 2}", 0);
    assert_eval_eq("VertexCoreness[PathGraph[4]]", "{1, 1, 1, 1}", 0);
    assert_eval_eq("VertexCoreness[StarGraph[5]]", "{1, 1, 1, 1, 1}", 0);
    assert_eval_eq("VertexCoreness[Graph[{1,2,3},{}]]", "{0, 0, 0}", 0);
    /* Triangle with a pendant: the triangle is 2-core, the pendant 1-core. */
    assert_eval_eq("VertexCoreness[Graph[{1,2,3,4},{1<->2,2<->3,3<->1,3<->4}]]", "{2, 2, 2, 1}", 0);
    /* Direction ignored; max coreness is the degeneracy. */
    assert_eval_eq("VertexCoreness[Graph[{1,2,3},{1->2,2->3,3->1}]]", "{2, 2, 2}", 0);
    assert_eval_eq("Max[VertexCoreness[CompleteGraph[4]]]", "3", 0);
    /* Coreness ≥ k agrees with membership in KCoreComponents[g, k]. */
    assert_eval_eq("Count[VertexCoreness[Graph[{1,2,3,4},{1<->2,2<->3,3<->1,3<->4}]], _?(#>=2&)]", "3", 0);
    assert_eval_eq("Head[VertexCoreness[5]]", "VertexCoreness", 0);
}

static void test_transitive_reduction(void) {
    /* Removes shortcut/transitive edges. */
    assert_eval_eq("EdgeList[TransitiveReductionGraph[Graph[{1,2,3},{1->2,2->3,1->3}]]]", "{1 -> 2, 2 -> 3}", 0);
    assert_eval_eq("EdgeCount[TransitiveReductionGraph[Graph[{1,2,3},{1->2,2->3}]]]", "2", 0);
    assert_eval_eq("EdgeCount[TransitiveReductionGraph[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4}]]]", "3", 0);
    assert_eval_eq("EdgeCount[TransitiveReductionGraph[Graph[{1,2,3},{1->2,1->3,2->3}]]]", "2", 0);
    assert_eval_eq("EdgeCount[TransitiveReductionGraph[Graph[{1,2,3,4},{1->2,1->3,2->4,3->4}]]]", "4", 0);
    assert_eval_eq("VertexCount[TransitiveReductionGraph[Graph[{1,2,3,4},{1->2,2->3,3->4,1->4}]]]", "4", 0);
    /* Reduction preserves reachability (same transitive closure). */
    assert_eval_eq("EdgeCount[TransitiveClosure[TransitiveReductionGraph[Graph[{1,2,3},{1->2,2->3,1->3}]]]] == EdgeCount[TransitiveClosure[Graph[{1,2,3},{1->2,2->3,1->3}]]]", "True", 0);
    assert_eval_eq("EdgeCount[TransitiveReductionGraph[Graph[{1,2,3},{}]]]", "0", 0);
    /* Cyclic / undirected inputs stay unevaluated. */
    assert_eval_eq("Head[TransitiveReductionGraph[Graph[{1,2,3},{1->2,2->3,3->1}]]]", "TransitiveReductionGraph", 0);
    assert_eval_eq("Head[TransitiveReductionGraph[CycleGraph[4]]]", "TransitiveReductionGraph", 0);
    assert_eval_eq("Head[TransitiveReductionGraph[5]]", "TransitiveReductionGraph", 0);
}

static void test_subgraph(void) {
    /* Induced subgraph: keep vertices and edges with both endpoints kept. */
    assert_eval_eq("CompleteGraphQ[Subgraph[CompleteGraph[4], {1,2,3}]]", "True", 0);
    assert_eval_eq("VertexCount[Subgraph[CompleteGraph[4], {1,2,3}]]", "3", 0);
    assert_eval_eq("EdgeCount[Subgraph[CompleteGraph[4], {1,2,3}]]", "3", 0);
    assert_eval_eq("EdgeList[Subgraph[CycleGraph[5], {1,2,3}]]", "{1 <-> 2, 2 <-> 3}", 0);
    assert_eval_eq("VertexCount[Subgraph[CompleteGraph[4], {}]]", "0", 0);
    assert_eval_eq("EdgeCount[Subgraph[CompleteGraph[4], {1}]]", "0", 0);
    /* Order preserved, non-vertices ignored, duplicates collapsed. */
    assert_eval_eq("VertexList[Subgraph[CompleteGraph[4], {3,1}]]", "{3, 1}", 0);
    assert_eval_eq("VertexList[Subgraph[CompleteGraph[3], {1,2,9}]]", "{1, 2}", 0);
    assert_eval_eq("VertexCount[Subgraph[CompleteGraph[3], {1,1,2}]]", "2", 0);
    assert_eval_eq("EdgeList[Subgraph[Graph[{1,2,3},{1->2,2->3,1->3}], {1,2}]]", "{1 -> 2}", 0);
    assert_eval_eq("Head[Subgraph[5, {1}]]", "Subgraph", 0);
}

static void test_vertex_delete(void) {
    assert_eval_eq("CompleteGraphQ[VertexDelete[CompleteGraph[4], 1]]", "True", 0);
    assert_eval_eq("VertexList[VertexDelete[CompleteGraph[4], 1]]", "{2, 3, 4}", 0);
    assert_eval_eq("EdgeCount[VertexDelete[CompleteGraph[4], {1,2}]]", "1", 0);
    assert_eval_eq("VertexCount[VertexDelete[CompleteGraph[4], {1,2}]]", "2", 0);
    assert_eval_eq("EdgeCount[VertexDelete[PathGraph[3], 2]]", "0", 0);
    assert_eval_eq("VertexList[VertexDelete[PathGraph[3], 2]]", "{1, 3}", 0);
    assert_eval_eq("EdgeCount[VertexDelete[CycleGraph[4], 1]]", "2", 0);
    assert_eval_eq("EdgeList[VertexDelete[Graph[{1,2,3},{1->2,2->3}], 1]]", "{2 -> 3}", 0);
    assert_eval_eq("VertexCount[VertexDelete[CompleteGraph[3], 9]]", "3", 0);   /* non-vertex: unchanged */
    assert_eval_eq("VertexCount[VertexDelete[CompleteGraph[3], {1,2,3}]]", "0", 0);
    assert_eval_eq("Head[VertexDelete[5, 1]]", "VertexDelete", 0);
}

static void test_edge_delete(void) {
    /* Delete via <-> sugar; all vertices kept; matching is symmetric. */
    assert_eval_eq("EdgeCount[EdgeDelete[CycleGraph[4], 1<->2]]", "3", 0);
    assert_eval_eq("VertexCount[EdgeDelete[CycleGraph[4], 1<->2]]", "4", 0);
    assert_eval_eq("EdgeCount[EdgeDelete[CycleGraph[4], 2<->1]]", "3", 0);
    assert_eval_eq("EdgeCount[EdgeDelete[CycleGraph[4], UndirectedEdge[1,2]]]", "3", 0);
    assert_eval_eq("EdgeList[EdgeDelete[CycleGraph[4], 1<->2]]", "{2 <-> 3, 3 <-> 4, 4 <-> 1}", 0);
    /* Delete a list of edges; nonexistent edges are ignored. */
    assert_eval_eq("EdgeCount[EdgeDelete[CompleteGraph[4], {1<->2,3<->4}]]", "4", 0);
    assert_eval_eq("EdgeCount[EdgeDelete[PathGraph[3], 1<->3]]", "2", 0);
    /* Directed: delete via ->, not symmetric. */
    assert_eval_eq("EdgeList[EdgeDelete[Graph[{1,2,3},{1->2,2->3}], 1->2]]", "{2 -> 3}", 0);
    assert_eval_eq("EdgeCount[EdgeDelete[Graph[{1,2},{1->2}], 2->1]]", "1", 0);
    assert_eval_eq("GraphQ[EdgeDelete[CompleteGraph[4], 1<->2]]", "True", 0);
    assert_eval_eq("Head[EdgeDelete[5, 1<->2]]", "EdgeDelete", 0);
}

static void test_edge_add(void) {
    /* Adding a chord closes a path into a triangle. */
    assert_eval_eq("EdgeCount[EdgeAdd[PathGraph[3], 1<->3]]", "3", 0);
    assert_eval_eq("CompleteGraphQ[EdgeAdd[PathGraph[3], 1<->3]]", "True", 0);
    /* Existing / symmetric-existing edges are not duplicated. */
    assert_eval_eq("EdgeCount[EdgeAdd[CycleGraph[4], 1<->2]]", "4", 0);
    assert_eval_eq("EdgeCount[EdgeAdd[CycleGraph[4], 2<->1]]", "4", 0);
    /* Missing endpoints become new vertices. */
    assert_eval_eq("VertexCount[EdgeAdd[Graph[{1,2},{1<->2}], 2<->3]]", "3", 0);
    assert_eval_eq("EdgeCount[EdgeAdd[Graph[{1,2},{1<->2}], 2<->3]]", "2", 0);
    /* Self-loops skipped; list add; directed via ->. */
    assert_eval_eq("EdgeCount[EdgeAdd[PathGraph[3], 1<->1]]", "2", 0);
    assert_eval_eq("EdgeCount[EdgeAdd[Graph[{1,2,3,4},{}], {1<->2,3<->4}]]", "2", 0);
    assert_eval_eq("MemberQ[EdgeList[EdgeAdd[Graph[{1,2},{1->2}], 2->1]], DirectedEdge[2,1]]", "True", 0);
    assert_eval_eq("EdgeCount[EdgeAdd[Graph[{1,2},{1->2}], 2->1]]", "2", 0);
    assert_eval_eq("GraphQ[EdgeAdd[PathGraph[3], 1<->3]]", "True", 0);
    assert_eval_eq("Head[EdgeAdd[5, 1<->2]]", "EdgeAdd", 0);
}

static void test_vertex_add(void) {
    assert_eval_eq("VertexCount[VertexAdd[CompleteGraph[3], 4]]", "4", 0);
    assert_eval_eq("EdgeCount[VertexAdd[CompleteGraph[3], 4]]", "3", 0);   /* isolated */
    assert_eval_eq("VertexList[VertexAdd[CompleteGraph[3], 4]]", "{1, 2, 3, 4}", 0);
    assert_eval_eq("VertexCount[VertexAdd[Graph[{1,2},{1<->2}], {3,4}]]", "4", 0);
    assert_eval_eq("VertexCount[VertexAdd[CompleteGraph[3], 2]]", "3", 0);  /* duplicate ignored */
    assert_eval_eq("VertexList[VertexAdd[Graph[{1},{}], x]]", "{1, x}", 0);
    assert_eval_eq("Last[DegreeCentrality[VertexAdd[CompleteGraph[3], 4]]]", "0", 0);
    assert_eval_eq("GraphQ[VertexAdd[CompleteGraph[3], 4]]", "True", 0);
    assert_eval_eq("VertexCount[VertexAdd[Graph[{},{}], {1,2}]]", "2", 0);
    assert_eval_eq("Head[VertexAdd[5, 1]]", "VertexAdd", 0);
}

static void test_neighborhood_graph(void) {
    assert_eval_eq("CompleteGraphQ[NeighborhoodGraph[CompleteGraph[4], 1]]", "True", 0);
    assert_eval_eq("Sort[VertexList[NeighborhoodGraph[CycleGraph[5], 1]]]", "{1, 2, 5}", 0);
    assert_eval_eq("EdgeCount[NeighborhoodGraph[CycleGraph[5], 1]]", "2", 0);
    assert_eval_eq("VertexCount[NeighborhoodGraph[PathGraph[5], 3]]", "3", 0);
    assert_eval_eq("EdgeCount[NeighborhoodGraph[PathGraph[5], 3]]", "2", 0);
    assert_eval_eq("MemberQ[VertexList[NeighborhoodGraph[CycleGraph[5], 1]], 1]", "True", 0);
    /* k controls the radius: 0 → just v, larger k widens. */
    assert_eval_eq("VertexCount[NeighborhoodGraph[CompleteGraph[4], 1, 0]]", "1", 0);
    assert_eval_eq("EdgeCount[NeighborhoodGraph[CompleteGraph[4], 1, 0]]", "0", 0);
    assert_eval_eq("VertexCount[NeighborhoodGraph[PathGraph[5], 3, 2]]", "5", 0);
    /* Bad arguments stay unevaluated. */
    assert_eval_eq("Head[NeighborhoodGraph[CycleGraph[3], 9]]", "NeighborhoodGraph", 0);
    assert_eval_eq("Head[NeighborhoodGraph[5, 1]]", "NeighborhoodGraph", 0);
}

static void test_graph_disjoint_union(void) {
    assert_eval_eq("VertexCount[GraphDisjointUnion[CompleteGraph[3], CompleteGraph[3]]]", "6", 0);
    assert_eval_eq("EdgeCount[GraphDisjointUnion[CompleteGraph[3], CompleteGraph[3]]]", "6", 0);
    assert_eval_eq("Length[ConnectedComponents[GraphDisjointUnion[CompleteGraph[3], CompleteGraph[3]]]]", "2", 0);
    assert_eval_eq("VertexCount[GraphDisjointUnion[PathGraph[2], PathGraph[2]]]", "4", 0);
    assert_eval_eq("EdgeCount[GraphDisjointUnion[PathGraph[2], PathGraph[2]]]", "2", 0);
    assert_eval_eq("VertexList[GraphDisjointUnion[Graph[{a,b},{a<->b}], Graph[{x},{}]]]", "{1, 2, 3}", 0);
    assert_eval_eq("ConnectedGraphQ[GraphDisjointUnion[CompleteGraph[3], CompleteGraph[3]]]", "False", 0);
    assert_eval_eq("EdgeCount[GraphDisjointUnion[CycleGraph[4], PathGraph[3]]] == 4+2", "True", 0);
    /* No cross edges (unlike GraphJoin); directed edges preserved. */
    assert_eval_eq("EdgeCount[GraphDisjointUnion[CompleteGraph[2], CompleteGraph[2]]] < EdgeCount[GraphJoin[CompleteGraph[2], CompleteGraph[2]]]", "True", 0);
    assert_eval_eq("DirectedGraphQ[GraphDisjointUnion[Graph[{1,2},{1->2}], Graph[{1,2},{1->2}]]]", "True", 0);
    assert_eval_eq("Head[GraphDisjointUnion[5, CycleGraph[3]]]", "GraphDisjointUnion", 0);
}

static void test_edge_contract(void) {
    /* Contracting a triangle edge leaves a single edge. */
    assert_eval_eq("EdgeList[EdgeContract[Graph[{1,2,3},{1<->2,2<->3,3<->1}], 1<->2]]", "{1 <-> 3}", 0);
    assert_eval_eq("VertexList[EdgeContract[Graph[{1,2,3},{1<->2,2<->3,3<->1}], 1<->2]]", "{1, 3}", 0);
    assert_eval_eq("EdgeCount[EdgeContract[PathGraph[3], 1<->2]]", "1", 0);
    /* List form; equals the corresponding VertexContract. */
    assert_eval_eq("VertexCount[EdgeContract[CompleteGraph[4], {1,2}]]", "3", 0);
    assert_eval_eq("EdgeContract[CompleteGraph[4], 1<->2] === VertexContract[CompleteGraph[4], {1,2}]", "True", 0);
    assert_eval_eq("VertexCount[EdgeContract[CompleteGraph[4], 1<->2]]", "3", 0);
    assert_eval_eq("EdgeCount[EdgeContract[CompleteGraph[4], 1<->2]]", "3", 0);
    assert_eval_eq("EdgeCount[EdgeContract[Graph[{1,2,3},{1->2,2->3}], 1->2]]", "1", 0);
    /* Bad edge specs stay unevaluated. */
    assert_eval_eq("Head[EdgeContract[CycleGraph[3], 1<->1]]", "EdgeContract", 0);
    assert_eval_eq("Head[EdgeContract[CycleGraph[3], 1<->9]]", "EdgeContract", 0);
    assert_eval_eq("Head[EdgeContract[5, 1<->2]]", "EdgeContract", 0);
}

static void test_find_matching(void) {
    assert_eval_eq("Length[FindIndependentEdgeSet[CompleteGraph[4]]]", "2", 0);
    assert_eval_eq("Length[FindIndependentEdgeSet[CompleteGraph[6]]]", "3", 0);   /* perfect */
    assert_eval_eq("Length[FindIndependentEdgeSet[PathGraph[4]]]", "2", 0);
    assert_eval_eq("Length[FindIndependentEdgeSet[CycleGraph[4]]]", "2", 0);
    assert_eval_eq("Length[FindIndependentEdgeSet[CycleGraph[5]]]", "2", 0);
    assert_eval_eq("Length[FindIndependentEdgeSet[StarGraph[5]]]", "1", 0);
    assert_eval_eq("Length[FindIndependentEdgeSet[CompleteGraph[3]]]", "1", 0);
    assert_eval_eq("FindIndependentEdgeSet[Graph[{1,2,3},{}]]", "{}", 0);
    /* Independence: the matching covers 2*size distinct vertices. */
    assert_eval_eq("With[{mm=FindIndependentEdgeSet[PathGraph[4]]}, Length[Union[Flatten[mm/.UndirectedEdge->List]]] == 2*Length[mm]]", "True", 0);
    assert_eval_eq("Head[FindIndependentEdgeSet[5]]", "FindIndependentEdgeSet", 0);
}

static void test_find_dominating_set(void) {
    assert_eval_eq("FindDominatingSet[StarGraph[5]]", "{1}", 0);
    assert_eval_eq("Length[FindDominatingSet[CompleteGraph[5]]]", "1", 0);
    assert_eval_eq("Length[FindDominatingSet[PathGraph[4]]]", "2", 0);
    assert_eval_eq("Length[FindDominatingSet[CycleGraph[4]]]", "2", 0);
    assert_eval_eq("Length[FindDominatingSet[CycleGraph[6]]]", "2", 0);
    assert_eval_eq("Length[FindDominatingSet[CycleGraph[7]]]", "3", 0);
    assert_eval_eq("Length[FindDominatingSet[Graph[{1,2,3},{}]]]", "3", 0);   /* edgeless */
    assert_eval_eq("FindDominatingSet[Graph[{1},{}]]", "{1}", 0);
    assert_eval_eq("Length[FindDominatingSet[Graph[{1,2,3},{1->2,2->3,3->1}]]]", "1", 0);
    assert_eval_eq("Head[FindDominatingSet[5]]", "FindDominatingSet", 0);
}

static void test_find_edge_cover(void) {
    assert_eval_eq("Length[FindEdgeCover[PathGraph[4]]]", "2", 0);
    assert_eval_eq("Length[FindEdgeCover[StarGraph[5]]]", "4", 0);   /* every leaf edge */
    assert_eval_eq("Length[FindEdgeCover[CompleteGraph[4]]]", "2", 0);
    assert_eval_eq("Length[FindEdgeCover[CompleteGraph[3]]]", "2", 0);
    assert_eval_eq("Length[FindEdgeCover[CycleGraph[6]]]", "3", 0);
    assert_eval_eq("Length[FindEdgeCover[Graph[{1,2},{1<->2}]]]", "1", 0);
    /* No cover when there is an isolated vertex. */
    assert_eval_eq("FindEdgeCover[Graph[{1,2},{}]]", "{}", 0);
    /* Gallai: |min edge cover| = n - |max matching|; cover touches every vertex. */
    assert_eval_eq("Length[FindEdgeCover[PathGraph[5]]] == 5 - Length[FindIndependentEdgeSet[PathGraph[5]]]", "True", 0);
    assert_eval_eq("With[{ec=FindEdgeCover[CycleGraph[5]]}, Length[Union[Flatten[ec/.UndirectedEdge->List]]] == 5]", "True", 0);
    assert_eval_eq("Head[FindEdgeCover[5]]", "FindEdgeCover", 0);
}

static void test_find_vertex_coloring(void) {
    assert_eval_eq("FindVertexColoring[CompleteGraph[3]]", "{1, 2, 3}", 0);
    assert_eval_eq("FindVertexColoring[CycleGraph[4]]", "{1, 2, 1, 2}", 0);
    assert_eval_eq("FindVertexColoring[Graph[{1,2,3},{}]]", "{1, 1, 1}", 0);
    assert_eval_eq("FindVertexColoring[Graph[{1},{}]]", "{1}", 0);
    assert_eval_eq("Length[FindVertexColoring[CycleGraph[6]]]", "6", 0);
    /* Number of colours used equals the chromatic number. */
    assert_eval_eq("Max[FindVertexColoring[CycleGraph[5]]] == ChromaticNumber[CycleGraph[5]]", "True", 0);
    assert_eval_eq("Max[FindVertexColoring[CompleteGraph[4]]] == 4", "True", 0);
    assert_eval_eq("Max[FindVertexColoring[PathGraph[5]]]", "2", 0);
    assert_eval_eq("Max[FindVertexColoring[StarGraph[5]]]", "2", 0);
    assert_eval_eq("Length[Union[FindVertexColoring[CompleteGraph[4]]]]", "4", 0);
    assert_eval_eq("Head[FindVertexColoring[5]]", "FindVertexColoring", 0);
}

static void test_graph_assortativity(void) {
    /* Stars are perfectly disassortative; a path P4 is -1/2. */
    assert_eval_eq("GraphAssortativity[StarGraph[4]]", "-1", 0);
    assert_eval_eq("GraphAssortativity[StarGraph[5]]", "-1", 0);
    assert_eval_eq("GraphAssortativity[PathGraph[4]]", "-1/2", 0);
    /* Regular / edgeless graphs have Indeterminate assortativity (zero variance). */
    assert_eval_eq("GraphAssortativity[CycleGraph[5]]", "Indeterminate", 0);
    assert_eval_eq("GraphAssortativity[CompleteGraph[4]]", "Indeterminate", 0);
    assert_eval_eq("GraphAssortativity[Graph[{1,2,3},{}]]", "Indeterminate", 0);
    assert_eval_eq("GraphAssortativity[Graph[{1,2,3},{1->2,2->3,3->1}]]", "Indeterminate", 0);
    assert_eval_eq("With[{r=GraphAssortativity[PathGraph[5]]}, -1<=r<=1]", "True", 0);
    assert_eval_eq("Head[GraphAssortativity[5]]", "GraphAssortativity", 0);
}

static void test_incidence_list(void) {
    assert_eval_eq("IncidenceList[CycleGraph[4], 1]", "{1 <-> 2, 4 <-> 1}", 0);
    assert_eval_eq("Length[IncidenceList[CycleGraph[4], 1]]", "2", 0);
    assert_eval_eq("Length[IncidenceList[StarGraph[5], 1]]", "4", 0);
    assert_eval_eq("Length[IncidenceList[StarGraph[5], 2]]", "1", 0);
    assert_eval_eq("IncidenceList[CycleGraph[3], 9]", "{}", 0);
    assert_eval_eq("IncidenceList[Graph[{1,2},{}], 1]", "{}", 0);
    /* Directed: includes both in- and out-edges at the vertex. */
    assert_eval_eq("Length[IncidenceList[Graph[{1,2,3},{1->2,2->3}], 2]]", "2", 0);
    assert_eval_eq("IncidenceList[Graph[{1,2,3},{1->2,2->3}], 1]", "{1 -> 2}", 0);
    /* Incidence count matches the (undirected) vertex degree. */
    assert_eval_eq("Length[IncidenceList[CompleteGraph[4], 1]] == VertexDegree[CompleteGraph[4], 1]", "True", 0);
    assert_eval_eq("Head[IncidenceList[5, 1]]", "IncidenceList", 0);
}

static void test_vertex_components(void) {
    /* Directed path 1->2->3. */
    assert_eval_eq("VertexOutComponent[Graph[{1,2,3},{1->2,2->3}], 1]", "{1, 2, 3}", 0);
    assert_eval_eq("VertexOutComponent[Graph[{1,2,3},{1->2,2->3}], 2]", "{2, 3}", 0);
    assert_eval_eq("VertexOutComponent[Graph[{1,2,3},{1->2,2->3}], 3]", "{3}", 0);
    assert_eval_eq("VertexInComponent[Graph[{1,2,3},{1->2,2->3}], 3]", "{1, 2, 3}", 0);
    assert_eval_eq("VertexInComponent[Graph[{1,2,3},{1->2,2->3}], 1]", "{1}", 0);
    assert_eval_eq("VertexInComponent[Graph[{1,2,3},{1->2,2->3}], 2]", "{1, 2}", 0);
    /* Undirected: both give v's connected component. */
    assert_eval_eq("VertexOutComponent[CycleGraph[4], 1]", "{1, 2, 3, 4}", 0);
    assert_eval_eq("VertexInComponent[CycleGraph[4], 1]", "{1, 2, 3, 4}", 0);
    assert_eval_eq("VertexOutComponent[Graph[{1,2,3},{1<->2}], 1]", "{1, 2}", 0);
    /* Directed cycle: every vertex reaches every other. */
    assert_eval_eq("VertexOutComponent[Graph[{1,2,3},{1->2,2->3,3->1}], 2]", "{1, 2, 3}", 0);
    assert_eval_eq("Head[VertexOutComponent[CycleGraph[3], 9]]", "VertexOutComponent", 0);
    assert_eval_eq("Head[VertexInComponent[5, 1]]", "VertexInComponent", 0);
}

static void test_graph_periphery(void) {
    /* Periphery = eccentricity-maximizing vertices (dual of GraphCenter). */
    assert_eval_eq("GraphPeriphery[PathGraph[5]]", "{1, 5}", 0);
    assert_eval_eq("GraphCenter[PathGraph[5]]", "{3}", 0);
    assert_eval_eq("GraphPeriphery[PathGraph[6]]", "{1, 6}", 0);
    assert_eval_eq("Length[GraphPeriphery[PathGraph[4]]]", "2", 0);
    assert_eval_eq("GraphPeriphery[StarGraph[5]]", "{2, 3, 4, 5}", 0);
    /* Vertex-transitive graphs: every vertex is peripheral. */
    assert_eval_eq("GraphPeriphery[CycleGraph[5]]", "{1, 2, 3, 4, 5}", 0);
    assert_eval_eq("GraphPeriphery[CompleteGraph[4]]", "{1, 2, 3, 4}", 0);
    assert_eval_eq("GraphPeriphery[Graph[{1},{}]]", "{1}", 0);
    /* Disconnected: all vertices have infinite eccentricity. */
    assert_eval_eq("GraphPeriphery[Graph[{1,2,3},{1<->2}]]", "{1, 2, 3}", 0);
    assert_eval_eq("Head[GraphPeriphery[5]]", "GraphPeriphery", 0);
}

static void test_antiprism_graph(void) {
    /* A3 is the octahedron: 6 vertices, 12 edges, 4-regular. */
    assert_eval_eq("VertexCount[AntiprismGraph[3]]", "6", 0);
    assert_eval_eq("EdgeCount[AntiprismGraph[3]]", "12", 0);
    assert_eval_eq("RegularGraphQ[AntiprismGraph[3]]", "True", 0);
    assert_eval_eq("EdgeCount[AntiprismGraph[3]] == EdgeCount[TuranGraph[6,3]]", "True", 0);
    assert_eval_eq("VertexCount[AntiprismGraph[4]]", "8", 0);
    assert_eval_eq("EdgeCount[AntiprismGraph[4]]", "16", 0);
    assert_eval_eq("EdgeCount[AntiprismGraph[5]] == 4*5", "True", 0);
    assert_eval_eq("RegularGraphQ[AntiprismGraph[5]]", "True", 0);
    assert_eval_eq("First[DegreeCentrality[AntiprismGraph[4]]]", "4", 0);
    assert_eval_eq("ConnectedGraphQ[AntiprismGraph[4]]", "True", 0);
    assert_eval_eq("Head[AntiprismGraph[2]]", "AntiprismGraph", 0);
}

static void test_prism_graph(void) {
    assert_eval_eq("VertexCount[PrismGraph[3]]", "6", 0);
    assert_eval_eq("EdgeCount[PrismGraph[3]]", "9", 0);
    assert_eval_eq("RegularGraphQ[PrismGraph[3]]", "True", 0);
    assert_eval_eq("EdgeCount[PrismGraph[4]]", "12", 0);          /* cube */
    assert_eval_eq("BipartiteGraphQ[PrismGraph[4]]", "True", 0);
    assert_eval_eq("EdgeCount[PrismGraph[5]] == 3*5", "True", 0);
    assert_eval_eq("First[DegreeCentrality[PrismGraph[5]]]", "3", 0);
    assert_eval_eq("ConnectedGraphQ[PrismGraph[6]]", "True", 0);
    /* Isomorphic to GeneralizedPetersenGraph[n,1] (same edge count). */
    assert_eval_eq("EdgeCount[PrismGraph[5]] == EdgeCount[GeneralizedPetersenGraph[5,1]]", "True", 0);
    assert_eval_eq("Head[PrismGraph[2]]", "PrismGraph", 0);
}

static void test_sunlet_graph(void) {
    assert_eval_eq("VertexCount[SunletGraph[3]]", "6", 0);
    assert_eval_eq("EdgeCount[SunletGraph[3]]", "6", 0);
    assert_eval_eq("VertexCount[SunletGraph[5]]", "10", 0);
    assert_eval_eq("EdgeCount[SunletGraph[5]] == 2*5", "True", 0);
    assert_eval_eq("RegularGraphQ[SunletGraph[4]]", "False", 0);
    assert_eval_eq("First[DegreeCentrality[SunletGraph[4]]]", "3", 0);   /* cycle vertex */
    assert_eval_eq("Last[DegreeCentrality[SunletGraph[4]]]", "1", 0);    /* pendant */
    assert_eval_eq("ConnectedGraphQ[SunletGraph[5]]", "True", 0);
    assert_eval_eq("BipartiteGraphQ[SunletGraph[4]]", "True", 0);
    assert_eval_eq("BipartiteGraphQ[SunletGraph[3]]", "False", 0);
    assert_eval_eq("Head[SunletGraph[2]]", "SunletGraph", 0);
}

static void test_helm_graph(void) {
    assert_eval_eq("VertexCount[HelmGraph[3]]", "7", 0);
    assert_eval_eq("EdgeCount[HelmGraph[3]]", "9", 0);
    assert_eval_eq("VertexCount[HelmGraph[5]]", "11", 0);
    assert_eval_eq("EdgeCount[HelmGraph[5]] == 3*5", "True", 0);
    assert_eval_eq("First[DegreeCentrality[HelmGraph[4]]]", "4", 0);       /* hub */
    assert_eval_eq("DegreeCentrality[HelmGraph[4]][[2]]", "4", 0);         /* rim */
    assert_eval_eq("Last[DegreeCentrality[HelmGraph[4]]]", "1", 0);        /* pendant */
    assert_eval_eq("ConnectedGraphQ[HelmGraph[5]]", "True", 0);
    assert_eval_eq("RegularGraphQ[HelmGraph[4]]", "False", 0);
    assert_eval_eq("BipartiteGraphQ[HelmGraph[4]]", "False", 0);
    assert_eval_eq("Head[HelmGraph[2]]", "HelmGraph", 0);
}

static void test_gear_graph(void) {
    assert_eval_eq("VertexCount[GearGraph[3]]", "7", 0);
    assert_eval_eq("EdgeCount[GearGraph[3]]", "9", 0);
    assert_eval_eq("VertexCount[GearGraph[5]]", "11", 0);
    assert_eval_eq("EdgeCount[GearGraph[5]] == 3*5", "True", 0);
    /* Gear graphs are bipartite (χ = 2). */
    assert_eval_eq("BipartiteGraphQ[GearGraph[4]]", "True", 0);
    assert_eval_eq("BipartiteGraphQ[GearGraph[3]]", "True", 0);
    assert_eval_eq("ChromaticNumber[GearGraph[4]]", "2", 0);
    assert_eval_eq("First[DegreeCentrality[GearGraph[4]]]", "4", 0);   /* hub degree n */
    assert_eval_eq("ConnectedGraphQ[GearGraph[5]]", "True", 0);
    assert_eval_eq("RegularGraphQ[GearGraph[4]]", "False", 0);
    assert_eval_eq("Head[GearGraph[2]]", "GearGraph", 0);
}

static void test_edge_betweenness(void) {
    assert_eval_eq("EdgeBetweennessCentrality[PathGraph[3]]", "{2, 2}", 0);
    assert_eval_eq("EdgeBetweennessCentrality[PathGraph[4]]", "{3, 4, 3}", 0);
    assert_eval_eq("EdgeBetweennessCentrality[StarGraph[4]]", "{3, 3, 3}", 0);
    assert_eval_eq("EdgeBetweennessCentrality[CycleGraph[4]]", "{2, 2, 2, 2}", 0);
    assert_eval_eq("EdgeBetweennessCentrality[CycleGraph[5]]", "{3, 3, 3, 3, 3}", 0);
    assert_eval_eq("EdgeBetweennessCentrality[CompleteGraph[4]]", "{1, 1, 1, 1, 1, 1}", 0);
    assert_eval_eq("EdgeBetweennessCentrality[Graph[{1,2},{1<->2}]]", "{1}", 0);
    assert_eval_eq("EdgeBetweennessCentrality[Graph[{1,2,3},{}]]", "{}", 0);
    assert_eval_eq("EdgeBetweennessCentrality[Graph[{1,2,3},{1->2,2->3}]]", "{2, 2}", 0);
    assert_eval_eq("Length[EdgeBetweennessCentrality[CompleteGraph[5]]]", "10", 0);
    assert_eval_eq("Head[EdgeBetweennessCentrality[5]]", "EdgeBetweennessCentrality", 0);
}

static void test_dodecahedral_graph(void) {
    assert_eval_eq("VertexCount[DodecahedralGraph[]]", "20", 0);
    assert_eval_eq("EdgeCount[DodecahedralGraph[]]", "30", 0);
    assert_eval_eq("RegularGraphQ[DodecahedralGraph[]]", "True", 0);
    assert_eval_eq("First[DegreeCentrality[DodecahedralGraph[]]]", "3", 0);
    assert_eval_eq("ConnectedGraphQ[DodecahedralGraph[]]", "True", 0);
    assert_eval_eq("BipartiteGraphQ[DodecahedralGraph[]]", "False", 0);
    assert_eval_eq("ChromaticNumber[DodecahedralGraph[]]", "3", 0);
    assert_eval_eq("HamiltonianGraphQ[DodecahedralGraph[]]", "True", 0);
    assert_eval_eq("EdgeCount[DodecahedralGraph[]] == EdgeCount[GeneralizedPetersenGraph[10,2]]", "True", 0);
    assert_eval_eq("Head[DodecahedralGraph[5]]", "DodecahedralGraph", 0);
}

static void test_icosahedral_graph(void) {
    assert_eval_eq("VertexCount[IcosahedralGraph[]]", "12", 0);
    assert_eval_eq("EdgeCount[IcosahedralGraph[]]", "30", 0);
    assert_eval_eq("RegularGraphQ[IcosahedralGraph[]]", "True", 0);
    assert_eval_eq("Union[DegreeCentrality[IcosahedralGraph[]]]", "{5}", 0);   /* 5-regular */
    assert_eval_eq("ConnectedGraphQ[IcosahedralGraph[]]", "True", 0);
    assert_eval_eq("BipartiteGraphQ[IcosahedralGraph[]]", "False", 0);
    assert_eval_eq("ChromaticNumber[IcosahedralGraph[]]", "4", 0);
    assert_eval_eq("HamiltonianGraphQ[IcosahedralGraph[]]", "True", 0);
    assert_eval_eq("Head[IcosahedralGraph[5]]", "IcosahedralGraph", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_edge_sugar_normalization);
    TEST(test_vertex_derivation);
    TEST(test_summary_printing);
    TEST(test_graphq);
    TEST(test_rejections);
    TEST(test_inputform_roundtrip);
    TEST(test_query_builtins);
    TEST(test_query_undirected);
    TEST(test_matrix_views);
    TEST(test_generators);
    TEST(test_random_graph);
    TEST(test_shortest_path);
    TEST(test_components);
    TEST(test_spanning_and_connectivity);
    TEST(test_graphplot);
    TEST(test_graphplot_options);
    TEST(test_highlight_graph);
    TEST(test_graph3d);
    TEST(test_bipartite);
    TEST(test_metrics);
    TEST(test_acyclic);
    TEST(test_complement);
    TEST(test_kirchhoff);
    TEST(test_edge_connectivity);
    TEST(test_line_graph);
    TEST(test_eulerian);
    TEST(test_star_wheel);
    TEST(test_complete_multipartite);
    TEST(test_grid_hypercube);
    TEST(test_closeness);
    TEST(test_transitive_closure);
    TEST(test_betweenness);
    TEST(test_find_eulerian);
    TEST(test_find_hamiltonian);
    TEST(test_graph_power);
    TEST(test_find_cycle);
    TEST(test_graph_distance_matrix);
    TEST(test_graph_density);
    TEST(test_degree_centrality);
    TEST(test_find_hamiltonian_path);
    TEST(test_kcore_components);
    TEST(test_local_clustering);
    TEST(test_global_clustering);
    TEST(test_mean_clustering);
    TEST(test_find_clique);
    TEST(test_find_independent);
    TEST(test_find_vertex_cover);
    TEST(test_graph_reciprocity);
    TEST(test_chromatic_polynomial);
    TEST(test_chromatic_number);
    TEST(test_degree_sequence);
    TEST(test_tree_graph_q);
    TEST(test_strongly_connected_q);
    TEST(test_hamiltonian_graph_q);
    TEST(test_regular_graph_q);
    TEST(test_complete_graph_q);
    TEST(test_graph_union);
    TEST(test_graph_intersection);
    TEST(test_graph_difference);
    TEST(test_graph_reverse);
    TEST(test_path_graph_q);
    TEST(test_vertex_contract);
    TEST(test_pagerank_centrality);
    TEST(test_katz_centrality);
    TEST(test_graph_join);
    TEST(test_index_graph);
    TEST(test_empty_and_mixed_q);
    TEST(test_graph_product);
    TEST(test_turan_graph);
    TEST(test_complete_kary_tree);
    TEST(test_circulant_graph);
    TEST(test_ladder_graph);
    TEST(test_cocktail_party_graph);
    TEST(test_kneser_graph);
    TEST(test_generalized_petersen_graph);
    TEST(test_friendship_graph);
    TEST(test_vertex_coreness);
    TEST(test_transitive_reduction);
    TEST(test_subgraph);
    TEST(test_vertex_delete);
    TEST(test_edge_delete);
    TEST(test_edge_add);
    TEST(test_vertex_add);
    TEST(test_neighborhood_graph);
    TEST(test_graph_disjoint_union);
    TEST(test_edge_contract);
    TEST(test_find_matching);
    TEST(test_find_dominating_set);
    TEST(test_find_edge_cover);
    TEST(test_find_vertex_coloring);
    TEST(test_graph_assortativity);
    TEST(test_incidence_list);
    TEST(test_vertex_components);
    TEST(test_graph_periphery);
    TEST(test_antiprism_graph);
    TEST(test_prism_graph);
    TEST(test_sunlet_graph);
    TEST(test_helm_graph);
    TEST(test_gear_graph);
    TEST(test_edge_betweenness);
    TEST(test_dodecahedral_graph);
    TEST(test_icosahedral_graph);

    printf("All graph tests passed!\n");
    return 0;
}
