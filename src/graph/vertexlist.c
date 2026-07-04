/* vertexlist.c - VertexList[g]: the graph's vertices, in canonical order.
 * Thin reader over the canonical form; NULL (unevaluated) if g is not a graph.
 * Memory (SPEC section 4): returns a fresh copy; the evaluator frees res. */

#include "graph.h"
#include "expr.h"

Expr* builtin_vertex_list(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    return expr_copy(g->data.function.args[0]);
}
