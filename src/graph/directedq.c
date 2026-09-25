/* directedq.c - DirectedGraphQ[g]: True iff g is a valid graph with at least
 * one edge, all of them DirectedEdge. False otherwise -- including for an
 * edgeless graph, which (as in the Wolfram Language) counts as undirected, so
 * DirectedGraphQ and UndirectedGraphQ are never both True.
 * Memory (SPEC section 4): returns a fresh symbol; the evaluator frees res. */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"

Expr* builtin_directed_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    long nd = graph_directed_edge_count(g);   /* O(1) via the validated-graph memo */
    int directed = nd > 0
        && (size_t)nd == g->data.function.args[1]->data.function.arg_count;
    return expr_new_symbol(directed ? SYM_True : SYM_False);
}
