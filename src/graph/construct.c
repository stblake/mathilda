/* construct.c - builtin_graph: normalize, derive, validate, canonicalize.
 *
 * Accepts:
 *   Graph[edges]                        -- vertices derived from the edges (directed default)
 *   Graph[verts, edges]                 -- explicit vertex list
 *   Graph[verts, edges, EdgeWeight -> {w1, ..., wm}]
 *                                        -- explicit vertex list + per-edge weights, matched
 *                                           to `edges` by position; wrong length is malformed
 *                                           (left unevaluated), same as any other rejection
 *                                           below. Weighted graphs require the explicit-vertex
 *                                           form -- Graph[edges, EdgeWeight -> {...}] is not
 *                                           accepted (deliberately out of scope; see the plan).
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

/* True iff `opt` is Rule[EdgeWeight, List[...]] -- shape only, length is
 * checked by the caller once the edge count is known. */
static int is_edge_weight_rule(const Expr* opt) {
    if (!opt || opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2)
        return 0;
    const char* h = fn_head(opt);
    if (h != SYM_Rule) return 0;
    const Expr* key = opt->data.function.args[0];
    return key && key->type == EXPR_SYMBOL && key->data.symbol.name == SYM_EdgeWeight
        && graph_is_list(opt->data.function.args[1]);
}

/* Assemble the canonical Graph from `res`, or NULL if the shape is wrong or the
 * result would be invalid. */
static Expr* try_build_canonical(Expr* res) {
    size_t argc = res->data.function.arg_count;
    const Expr* verts_in = NULL;
    const Expr* edges_in = NULL;
    const Expr* weight_opt = NULL;

    if (argc == 1) {
        edges_in = res->data.function.args[0];
    } else if (argc == 2) {
        verts_in = res->data.function.args[0];
        edges_in = res->data.function.args[1];
    } else if (argc == 3) {
        verts_in = res->data.function.args[0];
        edges_in = res->data.function.args[1];
        weight_opt = res->data.function.args[2];
        if (!is_edge_weight_rule(weight_opt)) return NULL;
    } else {
        return NULL;
    }
    if (!graph_is_list(edges_in)) return NULL;
    if (verts_in && !graph_is_list(verts_in)) return NULL;
    if (weight_opt && weight_opt->data.function.args[1]->data.function.arg_count
                       != edges_in->data.function.arg_count)
        return NULL;                          /* weight/edge count mismatch */

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
        /* At most 2 distinct new vertices per edge. First-appearance order is
         * preserved by appending on first sight; membership goes through a hash
         * index, because the linear expr_eq rescan it replaces made deriving the
         * vertex list O(E*V) -- seconds for a 40000-edge graph. */
        verts = (ne > 0) ? calloc(ne * 2, sizeof(Expr*)) : NULL;
        if (ne > 0 && !verts) goto fail_edges;
        GraphVIdx* seen = graph_vidx_new(ne * 2);
        if (!seen) { free(verts); goto fail_edges; }
        for (size_t i = 0; i < ne; i++) {
            for (int k = 0; k < 2; k++) {
                Expr* ep = edges[i]->data.function.args[k];
                /* Keys borrow `ep`, which lives in edges[i] until ownership
                 * moves into the result below -- outliving the index. */
                if (graph_vidx_put(seen, ep, (int)nv)) verts[nv++] = expr_copy(ep);
            }
        }
        graph_vidx_free(seen);
    }

    /* 3. Assemble candidate Graph[List verts, List edges(, EdgeWeight -> List w)]
     * (moves ownership). */
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, nv);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), edges, ne);
    free(verts);
    free(edges);
    Expr* g;
    if (weight_opt) {
        const Expr* win = weight_opt->data.function.args[1];
        size_t nw = win->data.function.arg_count;
        Expr** weights = (nw > 0) ? calloc(nw, sizeof(Expr*)) : NULL;
        if (nw > 0 && !weights) { expr_free(vlist); expr_free(elist); return NULL; }
        for (size_t i = 0; i < nw; i++)
            weights[i] = expr_copy(win->data.function.args[i]);
        Expr* wlist = expr_new_function(expr_new_symbol(SYM_List), weights, nw);
        free(weights);
        Expr* wargs[2] = { expr_new_symbol(SYM_EdgeWeight), wlist };
        Expr* wrule = expr_new_function(expr_new_symbol(SYM_Rule), wargs, 2);
        Expr* gargs[3] = { vlist, elist, wrule };
        g = expr_new_function(expr_new_symbol(SYM_Graph), gargs, 3);
    } else {
        Expr* gargs[2] = { vlist, elist };
        g = expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
    }

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
