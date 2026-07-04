/* vertexcontract.c - VertexContract[g, {v1, v2, ...}]: merge the listed vertices
 * into a single vertex (the first of the list), redirecting every incident edge
 * to the representative, deleting the resulting self-loops, and collapsing
 * parallel edges.
 *
 * The representative keeps its position; the other listed vertices are removed.
 * Each edge endpoint that names a merged vertex is rewritten to the
 * representative; an edge whose endpoints then coincide is a self-loop and is
 * dropped; remaining edges are deduplicated with edge-kind-aware,
 * symmetric-for-undirected equality. O((V+E)*k + E^2) at small-graph scale.
 * Contracting the two endpoints of an edge realises edge contraction (e.g. a
 * triangle becomes a single edge).
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL unless
 * arg1 is a valid graph and arg2 a List of its vertices.
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

/* Is v one of the merged vertices (list elements)? */
static int in_set(const Expr* set, const Expr* v) {
    for (size_t i = 0; i < set->data.function.arg_count; i++)
        if (expr_eq(set->data.function.args[i], v)) return 1;
    return 0;
}

Expr* builtin_vertex_contract(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    const Expr* set = res->data.function.args[1];
    if (!graph_is_valid(g) || !graph_is_list(set) || set->data.function.arg_count == 0)
        return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    size_t n = verts->data.function.arg_count;
    size_t m = edges->data.function.arg_count;
    const Expr* rep = set->data.function.args[0];

    /* Every listed vertex must be a vertex of g. */
    for (size_t i = 0; i < set->data.function.arg_count; i++)
        if (graph_vertex_index(verts, set->data.function.args[i]) < 0) return NULL;

    /* New vertices: keep rep and every non-merged vertex (order preserved). */
    Expr** vc = malloc((n > 0 ? n : 1) * sizeof(Expr*));
    size_t nv = 0;
    for (size_t i = 0; i < n; i++) {
        const Expr* v = verts->data.function.args[i];
        if (in_set(set, v) && !expr_eq(v, rep)) continue;   /* absorbed */
        vc[nv++] = expr_copy((Expr*)v);
    }

    /* Remap edges, drop self-loops, dedup. */
    Expr** ec = malloc((m > 0 ? m : 1) * sizeof(Expr*));
    size_t me = 0;
    for (size_t i = 0; i < m; i++) {
        const Expr* e = edges->data.function.args[i];
        const char* kind = graph_edge_kind(e);
        const Expr* a = e->data.function.args[0];
        const Expr* b = e->data.function.args[1];
        const Expr* na = in_set(set, a) ? rep : a;
        const Expr* nb = in_set(set, b) ? rep : b;
        if (expr_eq(na, nb)) continue;                       /* self-loop */
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
