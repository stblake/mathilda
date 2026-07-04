/* subgraph.c - Subgraph[g, {v1, v2, ...}]: the subgraph of g induced by the given
 * vertices -- those vertices (in the order listed, de-duplicated, restricted to
 * vertices actually in g) together with exactly the edges of g whose BOTH
 * endpoints are among them.
 *
 * O(V_sub + E) with linear membership tests (small-graph scale). Edge kinds are
 * preserved. Subgraph[K4, {1,2,3}] is K3; Subgraph[CycleGraph[5], {1,2,3}] is the
 * path 1-2-3; an empty vertex list gives the empty graph.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL unless
 * g is a valid graph and the second argument is a List.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_subgraph(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    const Expr* sel = res->data.function.args[1];
    if (!graph_is_valid(g) || !graph_is_list(sel)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t m = edges->data.function.arg_count;

    /* keep[i] = 1 if vertex i of g is selected. Also record selection order. */
    char* keep = (n > 0) ? calloc((size_t)n, 1) : NULL;
    Expr** vc = malloc((sel->data.function.arg_count > 0 ? sel->data.function.arg_count : 1) * sizeof(Expr*));
    size_t nv = 0;
    for (size_t i = 0; i < sel->data.function.arg_count; i++) {
        int vi = graph_vertex_index(verts, sel->data.function.args[i]);
        if (vi < 0 || keep[vi]) continue;          /* not in g, or duplicate */
        keep[vi] = 1;
        vc[nv++] = expr_copy(verts->data.function.args[vi]);
    }

    /* Keep edges with both endpoints selected. */
    Expr** ec = (m > 0) ? malloc(m * sizeof(Expr*)) : NULL;
    size_t me = 0;
    for (size_t k = 0; k < m; k++) {
        const Expr* e = edges->data.function.args[k];
        int a = graph_vertex_index(verts, e->data.function.args[0]);
        int b = graph_vertex_index(verts, e->data.function.args[1]);
        if (a >= 0 && b >= 0 && keep[a] && keep[b]) ec[me++] = expr_copy((Expr*)e);
    }
    free(keep);

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, nv);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
