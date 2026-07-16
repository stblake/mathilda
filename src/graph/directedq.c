/* directedq.c - DirectedGraphQ[g]: True iff g is a valid graph whose edges are
 * all DirectedEdge (vacuously True for an edgeless graph). False otherwise.
 * Memory (SPEC section 4): returns a fresh symbol; the evaluator frees res. */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"

Expr* builtin_directed_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return expr_new_symbol(SYM_False);

    const Expr* edges = g->data.function.args[1];
    for (size_t i = 0; i < edges->data.function.arg_count; i++) {
        if (graph_edge_kind(edges->data.function.args[i]) != SYM_DirectedEdge)
            return expr_new_symbol(SYM_False);
    }
    return expr_new_symbol(SYM_True);
}
