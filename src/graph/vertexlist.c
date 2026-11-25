/* vertexlist.c - VertexList[g] / VertexList[g, patt]: the graph's vertices, in
 * canonical order (those matching patt, for the two-argument form).
 * g may also be a list of edges ({1 -> 2, ...}, see graph_coerce). NULL
 * (unevaluated) if g is not a graph (a Hypergraph goes to its own reader).
 * Memory (SPEC section 4): returns a fresh copy; the evaluator frees res. */

#include "graph.h"
#include "expr.h"
#include "graph_hyper.h"

Expr* builtin_vertex_list(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    Expr* g = graph_coerce(res->data.function.args[0]);
    if (!g) return (argc == 1) ? hyp_vertex_list(res) : NULL;
    Expr* out = (argc == 2)
        ? graph_filter_pattern(g->data.function.args[0], res->data.function.args[1], 0, 0)
        : expr_copy(g->data.function.args[0]);
    expr_free(g);
    return out;
}
