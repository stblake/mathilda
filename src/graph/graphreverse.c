/* graphreverse.c - ReverseGraph[g]: g with the direction of every edge reversed.
 * A DirectedEdge[a, b] becomes DirectedEdge[b, a]; an UndirectedEdge is left
 * unchanged (it has no orientation). Vertices are preserved.
 *
 * The reverse (transpose) graph swaps successors and predecessors, so it turns
 * out-degree into in-degree, reachability into reverse reachability, and is an
 * involution: ReverseGraph[ReverseGraph[g]] === g. For an undirected graph it is
 * the identity. O(V + E), returns a canonical Graph.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL on a
 * non-graph argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_graph_reverse(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    size_t n = verts->data.function.arg_count;
    size_t m = edges->data.function.arg_count;

    Expr** vc = malloc((n > 0 ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) vc[i] = expr_copy(verts->data.function.args[i]);

    Expr** ec = malloc((m > 0 ? m : 1) * sizeof(Expr*));
    for (size_t i = 0; i < m; i++) {
        const Expr* e = edges->data.function.args[i];
        const char* kind = graph_edge_kind(e);
        const Expr* a = e->data.function.args[0];
        const Expr* b = e->data.function.args[1];
        if (kind == SYM_DirectedEdge) {          /* swap endpoints */
            Expr* args[2] = { expr_copy((Expr*)b), expr_copy((Expr*)a) };
            ec[i] = expr_new_function(expr_new_symbol(SYM_DirectedEdge), args, 2);
        } else {                                  /* undirected: unchanged */
            ec[i] = expr_copy((Expr*)e);
        }
    }

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, n);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, m);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
