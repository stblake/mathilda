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
#include "graph.h"
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

    /* --- RandomGraph[{n,m},k]: k independent graphs (RG-1) --- */
    assert_eval_eq("Length[RandomGraph[{6, 5}, 3]]", "3", 0);
    assert_eval_eq("And @@ (GraphQ /@ RandomGraph[{6, 5}, 3])", "True", 0);
    assert_eval_eq("Union[VertexCount /@ RandomGraph[{6, 5}, 3]]", "{6}", 0);
    assert_eval_eq("Union[EdgeCount /@ RandomGraph[{6, 5}, 3]]", "{5}", 0);
    /* k = 1 is a one-element list, not a bare graph. */
    assert_eval_eq("Length[RandomGraph[{6, 5}, 1]]", "1", 0);
    assert_eval_eq("RandomGraph[{6, 5}, 0]", "{}", 0);
    /* Bad k: silently unevaluated, per the src/random.c convention. */
    assert_eval_eq("Head[RandomGraph[{6, 5}, -1]]", "RandomGraph", 0);
    assert_eval_eq("Head[RandomGraph[{6, 5}, 2.5]]", "RandomGraph", 0);
    assert_eval_eq("Head[RandomGraph[{6, 5}, q]]", "RandomGraph", 0);
    /* m out of range fails identically at either arity. */
    assert_eval_eq("Head[RandomGraph[{3, 10}, 2]]", "RandomGraph", 0);
    /* Independence: 5 draws from C(28,4) are not all the same edge set. */
    assert_eval_eq("Length[Union[EdgeList /@ RandomGraph[{8, 4}, 5]]] > 1", "True", 0);
    /* n <= 1 and m = 0: edgeless graphs, no longer unevaluated. */
    assert_eval_eq("VertexCount[RandomGraph[{0, 0}]]", "0", 0);
    assert_eval_eq("EdgeCount[RandomGraph[{1, 0}]]", "0", 0);
    assert_eval_eq("Union[EdgeCount /@ RandomGraph[{5, 0}, 2]]", "{0}", 0);
    /* An absurd n or k is unevaluated, not a crash. n = 2^32 + 1 is the witness
     * that matters: n(n-1) wraps to 2^32 in 64 bits, so a guard applied to the
     * product instead of to n admits a small maxe and then overruns the
     * candidate buffer. n = 2^62 happens to wrap above the bound and so passes
     * either way. */
    assert_eval_eq("Head[RandomGraph[{4294967297, 1}, 1]]", "RandomGraph", 0);
    assert_eval_eq("Head[RandomGraph[{4294967297, 1}]]", "RandomGraph", 0);
    assert_eval_eq("Head[RandomGraph[{5, 2}, 4611686018427387904]]", "RandomGraph", 0);

    /* Determinism under a fixed seed, k form. */
    Expr* e2 = evaluate(parse_expression(
        "(SeedRandom[42]; EdgeList /@ RandomGraph[{6,5},3]) === "
        "(SeedRandom[42]; EdgeList /@ RandomGraph[{6,5},3])"));
    char* s2 = expr_to_string(e2);
    ASSERT(strcmp(s2, "True") == 0);
    free(s2);
    expr_free(e2);
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

/* ---- Vertex colouring internals (direct C, head not yet registered) -------
 * FindVertexColoring is deliberately unregistered until the search is proven
 * exact, so these cannot go through assert_eval_eq like every other test here.
 * Inputs are built the long way -- evaluate a generator expression, then
 * graph_build_adj -- which works because core_init() has run in main(). */

/* chi(expr) via the exact search, plus the backtracking-node count. */
static int coloured_chi(const char* src, long* steps_out, int** colour_out, int* n_out) {
    Expr* g = evaluate(parse_expression(src));
    GraphAdj* a = graph_build_adj(g);
    if (!a) { expr_free(g); return -1; }
    int* colour = calloc((size_t)(a->n > 0 ? a->n : 1), sizeof(int));
    int chi = fvc_search(a, colour, steps_out);
    int n = a->n;
    /* Verify properness here rather than in each caller: an exact search that
     * returns the right COUNT but an invalid assignment would otherwise pass. */
    for (int u = 0; u < n; u++) {
        for (int j = 0; j < a->outdeg[u]; j++)
            ASSERT_MSG(colour[u] != colour[a->out[u][j]], "adjacent vertices share a colour");
        ASSERT_MSG(colour[u] >= 1 && colour[u] <= chi, "colour outside 1..chi");
    }
    graph_adj_free(a);
    expr_free(g);
    if (colour_out) *colour_out = colour; else free(colour);
    if (n_out) *n_out = n;
    return chi;
}

static void test_vertex_coloring_internals(void) {
    long steps = 0;

    /* A complete graph needs one colour per vertex. */
    ASSERT_MSG(coloured_chi("CompleteGraph[5]", &steps, NULL, NULL) == 5,
                "chi(K5) should be 5");

    /* Even cycle is bipartite, odd cycle is not. */
    ASSERT_MSG(coloured_chi("CycleGraph[6]", &steps, NULL, NULL) == 2,
                "chi(C6) should be 2");
    ASSERT_MSG(coloured_chi("CycleGraph[5]", &steps, NULL, NULL) == 3,
                "chi(C5) should be 3");

    /* K_{2,2}: the discriminating case. A greedy pass on an unlucky vertex
     * order returns 3; only a minimal search returns 2. If minimality
     * regresses, this is the row that fails first. */
    ASSERT_MSG(coloured_chi("Graph[{1,2,3,4},{1<->3,1<->4,2<->3,2<->4}]",
                             &steps, NULL, NULL) == 2,
                "chi(K_{2,2}) should be 2, not a greedy 3");

    /* CompleteGraph[128] is UNDER the cap and therefore accepted. The clique
     * lower bound equals the DSATUR upper bound, so it must answer with ZERO
     * backtracking nodes -- asserting the counter, not merely that it is fast.
     * Without the lower bound this would refute k = 1..127 first, i.e. hang. */
    steps = 12345;
    ASSERT_MSG(coloured_chi("CompleteGraph[128]", &steps, NULL, NULL) == 128,
                "chi(K128) should be 128");
    ASSERT_MSG(steps == 0, "K128 must take zero search steps (lb == ub)");

    /* A sparse graph at the cap: bipartite, so also a bounds short-circuit. */
    ASSERT_MSG(coloured_chi("CycleGraph[128]", &steps, NULL, NULL) == 2,
                "chi(C128) should be 2");

    /* Degenerate shapes the builtin will hand straight to the search. */
    ASSERT_MSG(coloured_chi("Graph[{1},{}]", &steps, NULL, NULL) == 1,
                "chi(single vertex) should be 1");
    ASSERT_MSG(coloured_chi("Graph[{1,2,3,4},{}]", &steps, NULL, NULL) == 1,
                "chi(edgeless) should be 1");

    /* Direction is ignored for adjacency, and disconnected components are
     * minimised over the whole graph rather than per component. */
    ASSERT_MSG(coloured_chi("Graph[{1,2},{1->2}]", &steps, NULL, NULL) == 2,
                "a directed edge still constrains both endpoints");
    ASSERT_MSG(coloured_chi("Graph[{1,2,3,4},{1<->2,3<->4}]", &steps, NULL, NULL) == 2,
                "chi(two disjoint edges) should be 2");

    /* The bounds must bracket the true value on a case where they differ. */
    {
        Expr* g = evaluate(parse_expression("CycleGraph[5]"));
        GraphAdj* a = graph_build_adj(g);
        int* c = calloc(5, sizeof(int));
        int ub = fvc_dsatur_bound(a, c);
        int lb = fvc_clique_bound(a);
        ASSERT_MSG(lb == 2, "greedy clique in C5 should be an edge, size 2");
        ASSERT_MSG(ub >= 3, "DSATUR on C5 cannot beat chi = 3");
        free(c); graph_adj_free(a); expr_free(g);
    }

    /* Exported-symbol precondition (RG-2 adversarial finding #2). fvc_bb's
     * seen[]/forbid[] are stack arrays sized to FVC_MAX_VERTICES, and the cap is
     * enforced in the head -- not in fvc_search, which graph.h exports. A direct
     * caller passing n > cap must be refused (0), never overflow the buffers.
     * CompleteGraph[129] is one vertex over the cap and drives exactly that. */
    {
        Expr* g = evaluate(parse_expression("CompleteGraph[129]"));
        GraphAdj* a = graph_build_adj(g);
        ASSERT(a != NULL && a->n == 129);
        int* c = calloc((size_t)a->n, sizeof(int));
        steps = 999;
        ASSERT_MSG(fvc_search(a, c, &steps) == 0,
                    "fvc_search must refuse n > FVC_MAX_VERTICES, not overflow");
        ASSERT_MSG(steps == 0, "a cap-refused search reports zero steps");
        free(c); graph_adj_free(a); expr_free(g);
    }
}

/* ---- FindVertexColoring: the registered head (AC-1 .. AC-18) --------------
 * Sits alongside test_vertex_coloring_internals rather than replacing it: that
 * one reaches the search directly and can assert the node counter, which is not
 * observable from the language; this one pins the head's language-level
 * contract. The three long-running rows (AC-10b/10d/10f) are in
 * tests/test_graph_slow.c, which is excluded from this suite. */
static void test_vertex_coloring(void) {
    /* AC-1 .. AC-4: chromatic numbers of the classic shapes. */
    assert_eval_eq("Max[FindVertexColoring[CompleteGraph[5]]]", "5", 0);
    assert_eval_eq("Max[FindVertexColoring[CycleGraph[6]]]", "2", 0);
    assert_eval_eq("Max[FindVertexColoring[CycleGraph[5]]]", "3", 0);
    assert_eval_eq("Max[FindVertexColoring[PathGraph[4]]]", "2", 0);

    /* AC-5: no edges, so every vertex takes colour 1. */
    assert_eval_eq("Union[FindVertexColoring[Graph[{1,2,3,4},{}]]]", "{1}", 0);

    /* AC-6: exactly one colour per vertex. */
    assert_eval_eq("Length[FindVertexColoring[CycleGraph[7]]]", "7", 0);

    /* AC-7: properness, checked through the language. Indexes the colour vector
     * by vertex LABEL, so it is valid only for a graph labelled 1..n in order --
     * which CycleGraph is. */
    assert_eval_eq("Module[{g=CycleGraph[5],c},c=FindVertexColoring[g];"
                   "And@@(c[[#[[1]]]]=!=c[[#[[2]]]]&/@(List@@@EdgeList[g]))]",
                   "True", 0);

    /* AC-8: the VertexList-order contract, on the path c-a-b with non-integer
     * vertices. Position-sensitive by construction: in {c,a,b} order a minimal
     * colouring is {1,2,1}; under a sorted {a,b,c} order `a` is the degree-2
     * middle vertex and the colouring is {1,2,2}. So col[[1]] === col[[3]] is
     * True only for the VertexList order -- a mere "first two differ" check
     * holds under both and catches nothing. */
    assert_eval_eq("Module[{g=Graph[{c,a,b},{c<->a,a<->b}],col},"
                   "col=FindVertexColoring[g];"
                   "{Length[col],col[[1]]===col[[3]]&&col[[1]]=!=col[[2]]}]",
                   "{3, True}", 0);

    /* AC-9: above FVC_MAX_VERTICES the head refuses rather than answering. */
    assert_eval_eq("Head[FindVertexColoring[CompleteGraph[129]]]",
                   "FindVertexColoring", 0);

    /* AC-10: exactly at the cap. An even cycle gives ub=2 = lb, so this
     * exercises the cap BOUNDARY only -- it searches nothing. The genuinely
     * searched instance at the cap is AC-10b, in the slow target. */
    assert_eval_eq("Length[FindVertexColoring[CycleGraph[128]]]", "128", 0);

    /* AC-11, AC-12: the degenerate shapes. */
    assert_eval_eq("FindVertexColoring[Graph[{1},{}]]", "{1}", 0);
    assert_eval_eq("FindVertexColoring[Graph[{},{}]]", "{}", 0);

    /* AC-13: not a graph. AC-14: a graph that BUILDS but has an edge endpoint
     * absent from the vertex list, so it actually reaches the head (a
     * self-loop would be rejected by Graph itself and prove nothing). */
    assert_eval_eq("Head[FindVertexColoring[5]]", "FindVertexColoring", 0);
    assert_eval_eq("Head[FindVertexColoring[Graph[{1,2},{1<->3}]]]",
                   "FindVertexColoring", 0);

    /* AC-15: an edge constrains both endpoints whichever way it points. */
    assert_eval_eq("Max[FindVertexColoring[Graph[{1,2},{1->2}]]]", "2", 0);

    /* AC-16: minimal over the whole graph, not per component. */
    assert_eval_eq("Max[FindVertexColoring[Graph[{1,2,3,4},{1<->2,3<->4}]]]", "2", 0);

    /* AC-17: form 3 is a Non-goal, so a second argument must leave the
     * expression unevaluated rather than be silently ignored. */
    assert_eval_eq("Head[FindVertexColoring[CycleGraph[4], 3]]",
                   "FindVertexColoring", 0);

    /* AC-18: K_{2,2}. Greedy on an unlucky order says 3; minimal says 2. */
    assert_eval_eq("Max[FindVertexColoring[Graph[{1,2,3,4},{1<->3,1<->4,2<->3,2<->4}]]]",
                   "2", 0);

    /* No options are registered, deliberately -- Method would be an option
     * surface built to express a choice that does not yet exist. */
    assert_eval_eq("Options[FindVertexColoring]", "{}", 0);

    /* The search is deterministic, so repeated calls agree exactly. */
    assert_eval_eq("FindVertexColoring[CycleGraph[5]] === FindVertexColoring[CycleGraph[5]]",
                   "True", 0);
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
    TEST(test_vertex_coloring_internals);
    TEST(test_vertex_coloring);

    printf("All graph tests passed!\n");
    return 0;
}
