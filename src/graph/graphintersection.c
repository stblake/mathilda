/* graphintersection.c - GraphIntersection[g1, g2]: the graph whose vertices are
 * common to both graphs and whose edges are present in both, matched by
 * structural identity.
 *
 * Vertices are those of g1 that also appear in g2 (g1 order preserved). An edge
 * of g1 is kept iff an equal edge appears in g2, where equality respects edge
 * kind and is symmetric in the endpoints for undirected edges. A kept edge's
 * endpoints are necessarily common vertices, so the result is a valid canonical
 * Graph. O(V1*V2 + E1*E2) from the linear membership scans (small-graph scale).
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL unless
 * both arguments are valid graphs.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Structural edge equality, symmetric for undirected edges. */
static int same_edge(const Expr* a, const Expr* b) {
    const char* ka = graph_edge_kind(a);
    const char* kb = graph_edge_kind(b);
    if (!ka || ka != kb) return 0;
    const Expr* a0 = a->data.function.args[0];
    const Expr* a1 = a->data.function.args[1];
    const Expr* b0 = b->data.function.args[0];
    const Expr* b1 = b->data.function.args[1];
    if (expr_eq(a0, b0) && expr_eq(a1, b1)) return 1;
    if (ka == SYM_UndirectedEdge && expr_eq(a0, b1) && expr_eq(a1, b0)) return 1;
    return 0;
}

Expr* builtin_graph_intersection(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g1 = res->data.function.args[0];
    const Expr* g2 = res->data.function.args[1];
    if (!graph_is_valid(g1) || !graph_is_valid(g2)) return NULL;

    const Expr* v1 = g1->data.function.args[0];
    const Expr* e1 = g1->data.function.args[1];
    const Expr* v2 = g2->data.function.args[0];
    const Expr* e2 = g2->data.function.args[1];
    size_t n1 = v1->data.function.arg_count;
    size_t m1 = e1->data.function.arg_count, m2 = e2->data.function.arg_count;

    /* Vertices of g1 that also occur in g2 (g1 order). */
    Expr** vc = malloc((n1 > 0 ? n1 : 1) * sizeof(Expr*));
    size_t nv = 0;
    for (size_t i = 0; i < n1; i++) {
        const Expr* v = v1->data.function.args[i];
        if (graph_vertex_index(v2, v) >= 0) vc[nv++] = expr_copy((Expr*)v);
    }

    /* Edges of g1 that also occur in g2. */
    Expr** ec = malloc((m1 > 0 ? m1 : 1) * sizeof(Expr*));
    size_t me = 0;
    for (size_t i = 0; i < m1; i++) {
        const Expr* e = e1->data.function.args[i];
        int both = 0;
        for (size_t j = 0; j < m2; j++)
            if (same_edge(e, e2->data.function.args[j])) { both = 1; break; }
        if (both) ec[me++] = expr_copy((Expr*)e);
    }

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, nv);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
