/* The edge-list form of a graph argument: {1 -> 2, 2 <-> 3} in place of a
 * Graph, for every head that takes a graph as its first argument.
 *
 * The reference accepts a list of edges (or edge sugar: Rule, TwoWayRule)
 * wherever a graph is expected -- GraphDistance[{1 <-> 2, 2 <-> 3}, 1, 3] is
 * 2.  A handful of heads read their argument through graph_coerce and
 * accepted it (VertexList, EdgeList, VertexCount, EdgeCount,
 * ConnectedComponents, TopologicalSort); the rest required an evaluated Graph
 * and stayed unevaluated.
 *
 * Rather than repeat the coercion in ~70 builtins, the heads listed below are
 * routed through one shim (graph_edgearg_init, run after every graph module
 * has registered).  It calls the head's own builtin unchanged, and only when
 * argument 1 is a non-empty List of edge-likes does it first build the Graph
 * (graph_coerce) and call the builtin on a COPY of the call with the Graph in
 * place.  If that declines, the original call stays unevaluated as written --
 * the user sees their own edge list, not a Graph they never typed.  The cost
 * on the ordinary path is one pointer scan of this table.
 *
 * Predicates (*Q) are deliberately absent: they answer about the expression
 * they are given, and a rule list is not a graph.
 */

#include "graph.h"
#include "symtab.h"
#include "sym_names.h"
#include "eval.h"
#include <stdlib.h>

/* Heads whose first argument is a graph. */
static const char* const EDGEARG_HEADS[] = {
    /* queries and matrix views */
    "AdjacencyList", "VertexDegree", "VertexInDegree", "VertexOutDegree",
    "EdgeWeight", "AdjacencyMatrix", "WeightedAdjacencyMatrix",
    "IncidenceMatrix", "KirchhoffMatrix", "EdgeRules", "VertexIndex",
    "EdgeIndex",
    /* paths, distances, components, spanning trees, colouring */
    "FindShortestPath", "GraphDistance", "GraphDistanceMatrix",
    "VertexEccentricity", "GraphDiameter", "GraphRadius", "GraphCenter",
    "GraphPeriphery", "MeanGraphDistance", "GraphDensity",
    "WeaklyConnectedComponents", "StronglyConnectedComponents",
    "FindSpanningTree", "VertexConnectivity", "FindVertexColoring",
    /* centrality and clustering */
    "DegreeCentrality", "ClosenessCentrality", "EccentricityCentrality",
    "BetweennessCentrality", "EdgeBetweennessCentrality",
    "PageRankCentrality", "EigenvectorCentrality", "KatzCentrality",
    "HITSCentrality", "GraphTriangleCount", "LocalClusteringCoefficient",
    "GlobalClusteringCoefficient", "MeanClusteringCoefficient",
    /* editing and transforms */
    "VertexAdd", "VertexDelete", "EdgeAdd", "EdgeDelete", "Subgraph",
    "NeighborhoodGraph", "VertexReplace", "IndexGraph", "GraphComplement",
    "ReverseGraph", "UndirectedGraph", "DirectedGraph", "LineGraph",
    /* cycles, paths, flows, covers, cliques, isomorphism */
    "FindEulerianCycle", "FindCycle", "FindPath", "FindMaximumFlow",
    "FindMinimumCut", "FindEdgeCut", "FindVertexCut", "EdgeConnectivity",
    "FindIndependentEdgeSet", "FindEdgeCover", "FindVertexCover",
    "FindIndependentVertexSet", "FindClique", "FindKClique",
    "FindHamiltonianCycle", "FindHamiltonianPath", "FindGraphIsomorphism",
    "CanonicalGraph", "GraphAutomorphismGroup",
};
#define N_EDGEARG (sizeof(EDGEARG_HEADS) / sizeof(EDGEARG_HEADS[0]))

/* Of those, the heads that also take an adjacency matrix for the graph
 * (FindMinimumCut[{{0, 1}, {1, 0}}]): read as AdjacencyGraph[m]. */
static const char* const MATRIXARG_HEADS[] = { "FindMinimumCut" };
#define N_MATRIXARG (sizeof(MATRIXARG_HEADS) / sizeof(MATRIXARG_HEADS[0]))
static SymbolDef* matrixarg_def[N_MATRIXARG];

/* The head's own builtin, keyed by the head's SymbolDef (filled at init). */
static SymbolDef*  edgearg_def[N_EDGEARG];
static BuiltinFunc edgearg_fn[N_EDGEARG];

/* True iff e is a non-empty List whose every element is a two-argument edge
 * or edge sugar -- the only shape the shim rewrites (no evaluation). */
static int is_edge_list(const Expr* e) {
    if (!graph_is_list(e) || e->data.function.arg_count == 0) return 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        const Expr* x = e->data.function.args[i];
        if (x->type != EXPR_FUNCTION || x->data.function.arg_count != 2 ||
            x->data.function.head->type != EXPR_SYMBOL) return 0;
        const char* h = x->data.function.head->data.symbol.name;
        if (h != SYM_Rule && h != SYM_TwoWayRule &&
            h != SYM_DirectedEdge && h != SYM_UndirectedEdge) return 0;
    }
    return 1;
}

/* True iff e is a non-empty List of non-empty Lists (a matrix candidate for
 * AdjacencyGraph, which validates the rest). */
static int is_matrix_like(const Expr* e) {
    if (!graph_is_list(e) || e->data.function.arg_count == 0) return 0;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!graph_is_list(e->data.function.args[i])) return 0;
    return 1;
}

static int takes_matrix(const SymbolDef* d) {
    for (size_t i = 0; i < N_MATRIXARG; i++) if (matrixarg_def[i] == d) return 1;
    return 0;
}

/* AdjacencyGraph[m] evaluated, if it is a valid graph; NULL otherwise. */
static Expr* graph_from_matrix(const Expr* m) {
    Expr* a = expr_copy((Expr*)m);
    Expr* call = expr_new_function(expr_new_symbol(SYM_AdjacencyGraph), &a, 1);
    Expr* g = evaluate(call);
    expr_free(call);
    if (!graph_is_valid(g)) { expr_free(g); return NULL; }
    return g;
}

static Expr* builtin_graph_edgearg_shim(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.head->type != EXPR_SYMBOL)
        return NULL;
    SymbolDef* d = symtab_node_def(res->data.function.head);
    BuiltinFunc fn = NULL;
    for (size_t i = 0; i < N_EDGEARG; i++)
        if (edgearg_def[i] == d) { fn = edgearg_fn[i]; break; }
    if (!fn) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 0) return fn(res);
    const Expr* a0 = res->data.function.args[0];
    int matrix = !is_edge_list(a0);
    if (matrix && !(takes_matrix(d) && is_matrix_like(a0))) return fn(res);

    Expr* g = matrix ? graph_from_matrix(a0) : graph_coerce(a0);
    if (!g) return NULL;
    Expr** args = malloc(n * sizeof(Expr*));
    if (!args) { expr_free(g); return NULL; }
    args[0] = g;
    for (size_t i = 1; i < n; i++) args[i] = expr_copy(res->data.function.args[i]);
    Expr* call = expr_new_function(expr_copy(res->data.function.head), args, n);
    free(args);
    /* The builtin's answer (if any) is independent of `call`, which the
     * evaluator would have freed after a success; free it here either way. */
    Expr* out = fn(call);
    expr_free(call);
    return out;
}

void graph_edgearg_init(void) {
    for (size_t i = 0; i < N_EDGEARG; i++) {
        SymbolDef* d = symtab_get_def(EDGEARG_HEADS[i]);
        edgearg_def[i] = d;
        edgearg_fn[i] = d ? d->builtin_func : NULL;
        if (d && d->builtin_func) d->builtin_func = builtin_graph_edgearg_shim;
    }
    for (size_t i = 0; i < N_MATRIXARG; i++)
        matrixarg_def[i] = symtab_get_def(MATRIXARG_HEADS[i]);
}
