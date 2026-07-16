/* counts.c - VertexCount[g] and EdgeCount[g]: cardinalities as integers.
 * Thin readers over the canonical form; NULL (unevaluated) if g is not a graph.
 * Memory (SPEC section 4): returns a fresh integer; the evaluator frees res. */

#include "graph.h"
#include "expr.h"

Expr* builtin_vertex_count(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    return expr_new_integer((int64_t)g->data.function.args[0]->data.function.arg_count);
}

Expr* builtin_edge_count(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    return expr_new_integer((int64_t)g->data.function.args[1]->data.function.arg_count);
}
