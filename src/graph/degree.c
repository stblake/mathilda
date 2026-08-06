/* degree.c - VertexDegree, VertexInDegree, VertexOutDegree.
 *
 * Each accepts VertexDegree[g] (a list of degrees, one per vertex in canonical
 * order) or VertexDegree[g, v] (the degree of a single vertex).
 *
 * Conventions (documented in docs/spec/builtins/graphs.md):
 *   - A DirectedEdge[a,b] adds 1 to out(a) and 1 to in(b); total degree of a
 *     vertex is in + out, so for a purely directed graph total = in + out.
 *   - An UndirectedEdge[a,b] is incident to both a and b, adding 1 to each of
 *     their in-, out-, and total degrees (so in = out = total for a purely
 *     undirected graph).
 * There are no self-loops, so no endpoint is double-counted within one edge.
 *
 * Memory (SPEC section 4): returns freshly-allocated integers/lists; the
 * evaluator frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

enum { DEG_TOTAL, DEG_IN, DEG_OUT };

/* Degree of vertex `v` in graph `g` under the given mode. */
static int64_t degree_of(const Expr* g, const Expr* v, int mode) {
    const Expr* edges = g->data.function.args[1];
    int64_t in = 0, out = 0, tot = 0;
    for (size_t i = 0; i < edges->data.function.arg_count; i++) {
        const Expr* e = edges->data.function.args[i];
        const char* kind = graph_edge_kind(e);
        const Expr* a = e->data.function.args[0];
        const Expr* b = e->data.function.args[1];
        int va = expr_eq(a, v);
        int vb = expr_eq(b, v);
        if (kind == SYM_UndirectedEdge) {
            /* Incident once: contributes to in, out, and total alike. */
            if (va || vb) { in++; out++; tot++; }
        } else { /* directed: total counts each incident edge once */
            if (va) { out++; tot++; }
            if (vb) { in++;  tot++; }
        }
    }
    if (mode == DEG_IN)  return in;
    if (mode == DEG_OUT) return out;
    return tot;
}

static Expr* degree_dispatch(Expr* res, int mode) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    const Expr* verts = g->data.function.args[0];

    if (argc == 2) {
        const Expr* v = res->data.function.args[1];
        if (graph_vertex_index(verts, v) < 0) return NULL;   /* v not a vertex */
        return expr_new_integer(degree_of(g, v, mode));
    }

    /* All degrees at once: one pass over the edges, accumulating into per-vertex
     * counters, rather than one O(E) scan per vertex (which was O(V*E) -- ~5 s
     * for 20000 vertices and 40000 edges). The per-edge rules below are the same
     * ones degree_of applies, just pushed to the endpoints instead of pulled. */
    size_t nv = verts->data.function.arg_count;
    Expr** out = (nv > 0) ? calloc(nv, sizeof(Expr*)) : NULL;
    if (nv > 0 && !out) return NULL;

    int64_t* deg = (nv > 0) ? calloc(nv, sizeof(int64_t)) : NULL;
    GraphVIdx* ix = graph_vidx_new(nv);
    if ((nv > 0 && !deg) || !ix) { free(out); free(deg); graph_vidx_free(ix); return NULL; }
    for (size_t i = 0; i < nv; i++)
        graph_vidx_put(ix, verts->data.function.args[i], (int)i);

    const Expr* edges = g->data.function.args[1];
    for (size_t k = 0; k < edges->data.function.arg_count; k++) {
        const Expr* e = edges->data.function.args[k];
        const char* kind = graph_edge_kind(e);
        int ia = graph_vidx_get(ix, e->data.function.args[0]);
        int ib = graph_vidx_get(ix, e->data.function.args[1]);
        if (ia < 0 || ib < 0) continue;                 /* validated: cannot happen */
        if (kind == SYM_UndirectedEdge) {
            /* Incident to both ends, and in == out == total. No self-loops, so
             * the two ends are distinct and neither is double-counted. */
            deg[ia]++; deg[ib]++;
        } else if (mode == DEG_IN)  { deg[ib]++; }
        else if (mode == DEG_OUT)   { deg[ia]++; }
        else                        { deg[ia]++; deg[ib]++; }   /* total */
    }

    /* Read through the index so that a repeated vertex reports the degree of its
     * first occurrence -- what the per-vertex expr_eq scan did. */
    for (size_t i = 0; i < nv; i++) {
        int at = graph_vidx_get(ix, verts->data.function.args[i]);
        out[i] = expr_new_integer(deg[at >= 0 ? at : (int)i]);
    }
    graph_vidx_free(ix);
    free(deg);

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), out, nv);
    free(out);
    return list;
}

Expr* builtin_vertex_degree(Expr* res)     { return degree_dispatch(res, DEG_TOTAL); }
Expr* builtin_vertex_in_degree(Expr* res)  { return degree_dispatch(res, DEG_IN);    }
Expr* builtin_vertex_out_degree(Expr* res) { return degree_dispatch(res, DEG_OUT);   }
