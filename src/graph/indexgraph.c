/* indexgraph.c - IndexGraph[g] / IndexGraph[g, k]: g with its vertices renamed to
 * consecutive integers starting at 1 (or at k), in their current order, with
 * every edge remapped accordingly and edge kinds preserved.
 *
 * Useful for normalising arbitrary vertex labels (symbols, strings, nested
 * expressions) to a canonical integer indexing before matrix work or
 * comparison. O(V + E) with the linear vertex-index lookup. Returns a canonical
 * Graph.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL on a
 * non-graph or a non-integer start.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_index_graph(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    long start = 1;
    if (argc == 2) {
        const Expr* s = res->data.function.args[1];
        if (s->type != EXPR_INTEGER) return NULL;
        start = (long)s->data.integer;
    }

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    size_t n = verts->data.function.arg_count;
    size_t m = edges->data.function.arg_count;

    Expr** vc = malloc((n > 0 ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) vc[i] = expr_new_integer(start + (long)i);

    Expr** ec = malloc((m > 0 ? m : 1) * sizeof(Expr*));
    for (size_t i = 0; i < m; i++) {
        const Expr* e = edges->data.function.args[i];
        const char* kind = graph_edge_kind(e);
        int ia = graph_vertex_index(verts, e->data.function.args[0]);
        int ib = graph_vertex_index(verts, e->data.function.args[1]);
        Expr* args[2] = { expr_new_integer(start + ia), expr_new_integer(start + ib) };
        ec[i] = expr_new_function(expr_new_symbol(kind), args, 2);
    }

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, n);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, m);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
