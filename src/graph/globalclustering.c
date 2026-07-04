/* globalclustering.c - GlobalClusteringCoefficient[g]: the graph's transitivity,
 * the fraction of connected vertex triples that are closed into triangles:
 *
 *   C = 3 * (#triangles) / (#connected triples) = (sum_v L_v) / (sum_v C(d_v,2))
 *
 * where L_v is the number of edges among v's neighbours (so sum_v L_v counts
 * each triangle three times = 3T) and C(d_v,2) is the number of paths of length
 * two centred at v. Rewritten to stay integral, C = 2*sum_v L_v / sum_v
 * d_v(d_v-1). When there are no connected triples (no vertex of degree >= 2) the
 * coefficient is 0 by convention.
 *
 * Distinct from the mean of the local coefficients: transitivity weights each
 * vertex by how many triples it anchors. Edge direction is ignored (underlying
 * undirected graph). Symmetric boolean adjacency, O(V * d_max^2). Exact rational.
 *
 * Memory (SPEC section 4): returns a fresh number; frees res. NULL on a
 * non-graph argument.
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

Expr* builtin_global_clustering_coefficient(Expr* res) {
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

    long sumL = 0, sumDD = 0;   /* sum_v L_v ; sum_v d_v(d_v-1) */
    int* nbr = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;
    for (int v = 0; v < n; v++) {
        int d = 0;
        for (int u = 0; u < n; u++) if (adj[(size_t)v * n + u]) nbr[d++] = u;
        sumDD += (long)d * (d - 1);
        for (int i = 0; i < d; i++)
            for (int j = i + 1; j < d; j++)
                if (adj[(size_t)nbr[i] * n + nbr[j]]) sumL++;
    }
    free(nbr); free(adj);

    if (sumDD == 0) return expr_new_integer(0);   /* no connected triples */
    return ratio(2 * sumL, sumDD);
}
