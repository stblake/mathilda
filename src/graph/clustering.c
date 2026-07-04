/* clustering.c - LocalClusteringCoefficient[g]: for each vertex, the fraction of
 * its neighbour pairs that are themselves adjacent -- the classic measure of how
 * tightly a vertex's neighbourhood is knit.
 *
 *   C_v = (edges among neighbours of v) / C(deg(v), 2) = 2 L_v / (d_v (d_v - 1))
 *
 * and C_v = 0 when deg(v) < 2 (no pair to close). Edge direction is ignored (the
 * coefficient is computed on the underlying undirected graph); results are exact
 * rationals in vertex order. A symmetric boolean adjacency makes the
 * neighbour-pair test O(1), so the whole vector costs O(V * d_max^2) <= O(V^3).
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include "eval.h"
#include <stdlib.h>

/* Exact num/den as a reduced value (integer or Rational). den > 0 required. */
static Expr* ratio(long num, long den) {
    Expr* pa[2] = { expr_new_integer(den), expr_new_integer(-1) };
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
    Expr* ta[2] = { expr_new_integer(num), inv };
    return evaluate(expr_new_function(expr_new_symbol(SYM_Times), ta, 2));
}

Expr* builtin_local_clustering_coefficient(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    char* adj = (n > 0) ? calloc((size_t)n * (size_t)n, 1) : NULL;
    for (size_t e = 0; e < ne; e++) {
        const Expr* ed = edges->data.function.args[e];
        int ia = graph_vertex_index(verts, ed->data.function.args[0]);
        int ib = graph_vertex_index(verts, ed->data.function.args[1]);
        if (ia < 0 || ib < 0 || ia == ib) continue;
        adj[(size_t)ia * n + ib] = adj[(size_t)ib * n + ia] = 1;
    }

    Expr** items = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    int* nbr = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;
    for (int v = 0; v < n; v++) {
        int d = 0;
        for (int u = 0; u < n; u++) if (adj[(size_t)v * n + u]) nbr[d++] = u;
        if (d < 2) { items[v] = expr_new_integer(0); continue; }
        long links = 0;
        for (int i = 0; i < d; i++)
            for (int j = i + 1; j < d; j++)
                if (adj[(size_t)nbr[i] * n + nbr[j]]) links++;
        items[v] = ratio(2 * links, (long)d * (d - 1));
    }
    free(nbr); free(adj);

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
    free(items);
    return out;
}
