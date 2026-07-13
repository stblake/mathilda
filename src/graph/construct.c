/* construct.c - builtin_graph: normalize, derive, validate, canonicalize.
 *
 * Accepts:
 *   Graph[edges]          -- vertices derived from the edges (directed default)
 *   Graph[verts, edges]   -- explicit vertex list
 *
 * Edge sugar is normalized on construction:
 *   Rule[u,v]        / u -> v    ->  DirectedEdge[u, v]
 *   TwoWayRule[u,v]  / u <-> v   ->  UndirectedEdge[u, v]
 *   DirectedEdge[u,v] / UndirectedEdge[u,v]   pass through unchanged
 *
 * The result is the canonical Graph[List[verts], List[edges]] with vertices in
 * first-appearance order (when derived). Malformed input -- 3-arg edges,
 * self-loops, parallel edges, or an edge endpoint absent from an explicit
 * vertex list -- leaves Graph[...] unevaluated (returns NULL).
 *
 * Memory (SPEC section 4): the canonical tree is built entirely from expr_copy
 * of the argument's parts, so `res` is never cannibalized. On success the
 * evaluator frees `res`; on NULL it retains it. The "already canonical" case
 * returns NULL so evaluation reaches a fixed point.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* head symbol of a function node, or NULL. */
static const char* fn_head(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || !e->data.function.head
        || e->data.function.head->type != EXPR_SYMBOL)
        return NULL;
    return e->data.function.head->data.symbol.name;
}

/* Build a fresh normalized edge (with copied endpoints) from any accepted edge
 * form, or NULL if `e` is not a recognizable 2-argument edge. Does not reject
 * self-loops -- the caller does, via graph_is_valid. */
static Expr* normalize_edge(const Expr* e) {
    const char* h = fn_head(e);
    if (!h || e->data.function.arg_count != 2) return NULL;

    const char* out_head;
    if (h == SYM_Rule || h == SYM_DirectedEdge)              out_head = SYM_DirectedEdge;
    else if (h == SYM_TwoWayRule || h == SYM_UndirectedEdge) out_head = SYM_UndirectedEdge;
    else return NULL;

    Expr* args[2];
    args[0] = expr_copy(e->data.function.args[0]);
    args[1] = expr_copy(e->data.function.args[1]);
    return expr_new_function(expr_new_symbol(out_head), args, 2);
}

/* Assemble the canonical Graph from `res`, or NULL if the shape is wrong or the
 * result would be invalid. */
static Expr* try_build_canonical(Expr* res) {
    size_t argc = res->data.function.arg_count;
    const Expr* verts_in = NULL;
    const Expr* edges_in = NULL;

    if (argc == 1) {
        edges_in = res->data.function.args[0];
    } else if (argc == 2) {
        verts_in = res->data.function.args[0];
        edges_in = res->data.function.args[1];
    } else {
        return NULL;
    }
    if (!graph_is_list(edges_in)) return NULL;
    if (verts_in && !graph_is_list(verts_in)) return NULL;

    size_t ne = edges_in->data.function.arg_count;

    /* 1. Normalize every edge. */
    Expr** edges = (ne > 0) ? calloc(ne, sizeof(Expr*)) : NULL;
    if (ne > 0 && !edges) return NULL;
    for (size_t i = 0; i < ne; i++) {
        edges[i] = normalize_edge(edges_in->data.function.args[i]);
        if (!edges[i]) {
            for (size_t j = 0; j < i; j++) expr_free(edges[j]);
            free(edges);
            return NULL;
        }
    }

    /* 2. Build the vertex list: copy the explicit one, or derive it from the
     *    normalized edges in first-appearance order. */
    Expr** verts = NULL;
    size_t nv = 0;
    if (verts_in) {
        nv = verts_in->data.function.arg_count;
        verts = (nv > 0) ? calloc(nv, sizeof(Expr*)) : NULL;
        if (nv > 0 && !verts) goto fail_edges;
        for (size_t i = 0; i < nv; i++)
            verts[i] = expr_copy(verts_in->data.function.args[i]);
    } else {
        /* At most 2 distinct new vertices per edge. */
        verts = (ne > 0) ? calloc(ne * 2, sizeof(Expr*)) : NULL;
        if (ne > 0 && !verts) goto fail_edges;
        for (size_t i = 0; i < ne; i++) {
            for (int k = 0; k < 2; k++) {
                Expr* ep = edges[i]->data.function.args[k];
                int seen = 0;
                for (size_t j = 0; j < nv; j++)
                    if (expr_eq(verts[j], ep)) { seen = 1; break; }
                if (!seen) verts[nv++] = expr_copy(ep);
            }
        }
    }

    /* 3. Assemble candidate Graph[List verts, List edges] (moves ownership). */
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, nv);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), edges, ne);
    free(verts);
    free(edges);
    Expr* gargs[2] = { vlist, elist };
    Expr* g = expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);

    /* 4. Validate (self-loops, parallel edges, endpoint membership). */
    if (!graph_is_valid(g)) { expr_free(g); return NULL; }
    return g;

fail_edges:
    for (size_t i = 0; i < ne; i++) expr_free(edges[i]);
    free(edges);
    return NULL;
}

Expr* builtin_graph(Expr* res) {
    Expr* canonical = try_build_canonical(res);
    if (!canonical) return NULL;                 /* malformed: leave unevaluated */
    if (expr_eq(canonical, res)) {               /* already canonical: fixed point */
        expr_free(canonical);
        return NULL;
    }
    return canonical;
}
