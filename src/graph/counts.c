/* counts.c - VertexCount[g] / VertexCount[g, patt] and EdgeCount[g] /
 * EdgeCount[g, patt]: cardinalities as integers.
 * Thin readers over the canonical form. g may also be a list of edges
 * ({1 -> 2, ...}, see graph_coerce). The pattern form counts the vertices or
 * edges matching patt (graph_filter_pattern). NULL (unevaluated) if g is not a
 * graph (a Hypergraph is handed to its own reader).
 * Memory (SPEC section 4): returns a fresh integer; the evaluator frees res. */

#include "graph.h"
#include "expr.h"
#include "graph_hyper.h"

/* Shared driver: which = 0 for vertices, 1 for edges. */
static Expr* count_driver(Expr* res, int which) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    Expr* g = graph_coerce(res->data.function.args[0]);
    if (!g) {
        if (argc != 1) return NULL;
        return which ? hyp_edge_count(res) : hyp_vertex_count(res);
    }
    const Expr* items = g->data.function.args[which];
    Expr* out = (argc == 2)
        ? graph_filter_pattern(items, res->data.function.args[1], which, 1)
        : expr_new_integer((int64_t)items->data.function.arg_count);
    expr_free(g);
    return out;
}

Expr* builtin_vertex_count(Expr* res) { return count_driver(res, 0); }
Expr* builtin_edge_count(Expr* res)   { return count_driver(res, 1); }
