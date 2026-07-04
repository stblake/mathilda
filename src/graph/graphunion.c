/* graphunion.c - GraphUnion[g1, g2]: the graph whose vertex set is the union of
 * the two vertex sets and whose edge set is the union of the two edge sets,
 * matched by vertex identity (structural equality). Duplicate vertices and
 * edges are collapsed.
 *
 * Vertices from g1 keep their order; new vertices from g2 are appended. Edges
 * likewise, with an edge from g2 added only if no equal edge is already present
 * -- where "equal" respects edge kind and, for undirected edges, is symmetric in
 * the endpoints. The result is returned as a canonical Graph[List, List].
 * O((V1+V2)^2 + (E1+E2)^2) from the linear membership scans (small-graph scale).
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

Expr* builtin_graph_union(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g1 = res->data.function.args[0];
    const Expr* g2 = res->data.function.args[1];
    if (!graph_is_valid(g1) || !graph_is_valid(g2)) return NULL;

    const Expr* v1 = g1->data.function.args[0];
    const Expr* e1 = g1->data.function.args[1];
    const Expr* v2 = g2->data.function.args[0];
    const Expr* e2 = g2->data.function.args[1];
    size_t n1 = v1->data.function.arg_count, n2 = v2->data.function.arg_count;
    size_t m1 = e1->data.function.arg_count, m2 = e2->data.function.arg_count;

    /* Union of vertices (g1 order, then new ones from g2). */
    const Expr** vp = malloc((n1 + n2) * sizeof(Expr*));
    size_t nv = 0;
    for (size_t i = 0; i < n1; i++) vp[nv++] = v1->data.function.args[i];
    for (size_t i = 0; i < n2; i++) {
        const Expr* v = v2->data.function.args[i];
        int seen = 0;
        for (size_t j = 0; j < nv; j++) if (expr_eq(vp[j], v)) { seen = 1; break; }
        if (!seen) vp[nv++] = v;
    }

    /* Union of edges. */
    const Expr** ep = malloc(((m1 + m2) > 0 ? (m1 + m2) : 1) * sizeof(Expr*));
    size_t me = 0;
    for (size_t i = 0; i < m1; i++) ep[me++] = e1->data.function.args[i];
    for (size_t i = 0; i < m2; i++) {
        const Expr* e = e2->data.function.args[i];
        int seen = 0;
        for (size_t j = 0; j < me; j++) if (same_edge(ep[j], e)) { seen = 1; break; }
        if (!seen) ep[me++] = e;
    }

    Expr** vc = malloc((nv > 0 ? nv : 1) * sizeof(Expr*));
    for (size_t i = 0; i < nv; i++) vc[i] = expr_copy((Expr*)vp[i]);
    Expr** ec = malloc((me > 0 ? me : 1) * sizeof(Expr*));
    for (size_t i = 0; i < me; i++) ec[i] = expr_copy((Expr*)ep[i]);
    free(vp); free(ep);

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, nv);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
