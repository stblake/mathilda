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

    printf("All graph tests passed!\n");
    return 0;
}
