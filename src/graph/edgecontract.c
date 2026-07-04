/* edgecontract.c - EdgeContract[g, e]: contract the edge e by merging its two
 * endpoints into a single vertex (the first endpoint), redirecting every incident
 * edge to it, deleting the resulting self-loop, and collapsing parallel edges.
 *
 * The edge may be given as DirectedEdge/UndirectedEdge, the sugar a->b / a<->b,
 * or a two-element list {a, b}; both endpoints must be vertices of g. This is the
 * edge-oriented view of vertex contraction: EdgeContract[g, {u,v}] equals
 * VertexContract[g, {u,v}]. O((V+E) + E^2) at small-graph scale; edge kinds are
 * preserved.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL unless
 * g is a valid graph and e names two vertices of g.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

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

/* Extract the two endpoints of an edge spec; returns 0 on failure. */
static int spec_endpoints(const Expr* e, const Expr** a, const Expr** b) {
    if (e->type != EXPR_FUNCTION || e->data.function.arg_count != 2 ||
        e->data.function.head->type != EXPR_SYMBOL) return 0;
    const char* h = e->data.function.head->data.symbol;
    if (h == SYM_DirectedEdge || h == SYM_UndirectedEdge ||
        h == SYM_Rule || h == SYM_TwoWayRule || h == SYM_List) {
        *a = e->data.function.args[0];
        *b = e->data.function.args[1];
        return 1;
    }
    return 0;
}

Expr* builtin_edge_contract(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    const Expr *ea, *eb;
    if (!spec_endpoints(res->data.function.args[1], &ea, &eb)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    size_t n = verts->data.function.arg_count;
    size_t m = edges->data.function.arg_count;
    if (graph_vertex_index(verts, ea) < 0 || graph_vertex_index(verts, eb) < 0) return NULL;
    if (expr_eq(ea, eb)) return NULL;                    /* not an edge */

    /* New vertices: keep everything except eb (absorbed into ea). */
    Expr** vc = malloc((n > 0 ? n : 1) * sizeof(Expr*));
    size_t nv = 0;
    for (size_t i = 0; i < n; i++) {
        const Expr* v = verts->data.function.args[i];
        if (expr_eq(v, eb)) continue;
        vc[nv++] = expr_copy((Expr*)v);
    }

    /* Remap edges (eb -> ea), drop self-loops, dedup. */
    Expr** ec = (m > 0) ? malloc(m * sizeof(Expr*)) : NULL;
    size_t me = 0;
    for (size_t k = 0; k < m; k++) {
        const Expr* e = edges->data.function.args[k];
        const char* kind = graph_edge_kind(e);
        const Expr* a = e->data.function.args[0];
        const Expr* b = e->data.function.args[1];
        const Expr* na = expr_eq(a, eb) ? ea : a;
        const Expr* nb = expr_eq(b, eb) ? ea : b;
        if (expr_eq(na, nb)) continue;                   /* self-loop */
        Expr* ne = expr_new_function(expr_new_symbol(kind),
                       (Expr*[]){ expr_copy((Expr*)na), expr_copy((Expr*)nb) }, 2);
        int dup = 0;
        for (size_t j = 0; j < me; j++) if (same_edge(ec[j], ne)) { dup = 1; break; }
        if (dup) expr_free(ne); else ec[me++] = ne;
    }

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, nv);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
