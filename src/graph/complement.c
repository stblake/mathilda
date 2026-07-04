/* complement.c - GraphComplement[g]: the graph on the same vertices whose edges
 * are exactly the non-edges of g.
 *
 *   undirected g -> undirected complement: every pair {i,j} (i<j) not already
 *                   joined becomes an UndirectedEdge (edgeless -> complete K_n).
 *   directed g   -> directed complement: every ordered pair (i,j), i != j, with
 *                   no i->j edge becomes a DirectedEdge.
 * No self-loops in either case. A graph with mixed edge kinds is complemented as
 * undirected (its underlying undirected structure).
 *
 * Uses an n x n membership matrix, so it is O(V^2) — the natural cost, since the
 * complement itself is dense. Returns the canonical Graph value directly.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph tree; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_graph_complement(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    /* Directed iff every edge is a DirectedEdge (edgeless is treated as
     * undirected, so its complement is the complete graph). */
    int directed = (ne > 0);
    for (size_t i = 0; i < ne; i++)
        if (graph_edge_kind(edges->data.function.args[i]) != SYM_DirectedEdge) { directed = 0; break; }

    char* present = (n > 0) ? calloc((size_t)n * (size_t)n, 1) : NULL;
    if (n > 0 && !present) return NULL;
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        int ia = graph_vertex_index(verts, e->data.function.args[0]);
        int ib = graph_vertex_index(verts, e->data.function.args[1]);
        if (ia < 0 || ib < 0) continue;
        present[(size_t)ia * n + ib] = 1;
        if (!directed) present[(size_t)ib * n + ia] = 1;
    }

    /* Upper bound on complement edges: n(n-1) directed, n(n-1)/2 undirected. */
    size_t cap = directed ? (size_t)n * (n > 0 ? n - 1 : 0)
                          : (size_t)n * (n > 0 ? n - 1 : 0) / 2;
    Expr** cedges = (cap > 0) ? calloc(cap, sizeof(Expr*)) : NULL;
    size_t m = 0;
    const char* ekind = directed ? SYM_DirectedEdge : SYM_UndirectedEdge;
    for (int i = 0; i < n; i++) {
        for (int j = directed ? 0 : i + 1; j < n; j++) {
            if (i == j) continue;
            if (present[(size_t)i * n + j]) continue;
            Expr* args[2] = { expr_copy(verts->data.function.args[i]),
                              expr_copy(verts->data.function.args[j]) };
            cedges[m++] = expr_new_function(expr_new_symbol(ekind), args, 2);
        }
    }
    free(present);

    Expr** vcopy = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    for (int i = 0; i < n; i++) vcopy[i] = expr_copy(verts->data.function.args[i]);
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vcopy, (size_t)n);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), cedges, m);
    free(vcopy); free(cedges);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
