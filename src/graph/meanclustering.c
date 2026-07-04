/* meanclustering.c - MeanClusteringCoefficient[g]: the average of the local
 * clustering coefficients, (1/n) * sum_v C_v, where C_v is the fraction of v's
 * neighbour pairs that are adjacent (0 for deg(v) < 2).
 *
 * This is Mean[LocalClusteringCoefficient[g]] -- every vertex counts equally,
 * including low-degree ones contributing 0. It differs from
 * GlobalClusteringCoefficient (transitivity), which weights vertices by the
 * number of triples they anchor. Edge direction is ignored; the value is an
 * exact rational (the per-vertex terms are summed and divided by n through the
 * evaluator, so the result is fully reduced). O(V * d_max^2).
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

Expr* builtin_mean_clustering_coefficient(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;
    if (n == 0) return expr_new_integer(0);   /* mean over no vertices */

    char* adj = calloc((size_t)n * (size_t)n, 1);
    for (size_t e = 0; e < ne; e++) {
        const Expr* ed = edges->data.function.args[e];
        int ia = graph_vertex_index(verts, ed->data.function.args[0]);
        int ib = graph_vertex_index(verts, ed->data.function.args[1]);
        if (ia < 0 || ib < 0 || ia == ib) continue;
        adj[(size_t)ia * n + ib] = adj[(size_t)ib * n + ia] = 1;
    }

    Expr** terms = calloc((size_t)n, sizeof(Expr*));
    int* nbr = malloc((size_t)n * sizeof(int));
    for (int v = 0; v < n; v++) {
        int d = 0;
        for (int u = 0; u < n; u++) if (adj[(size_t)v * n + u]) nbr[d++] = u;
        if (d < 2) { terms[v] = expr_new_integer(0); continue; }
        long links = 0;
        for (int i = 0; i < d; i++)
            for (int j = i + 1; j < d; j++)
                if (adj[(size_t)nbr[i] * n + nbr[j]]) links++;
        terms[v] = ratio(2 * links, (long)d * (d - 1));
    }
    free(nbr); free(adj);

    /* mean = (sum_v C_v) / n, reduced by the evaluator. */
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, (size_t)n);
    free(terms);
    Expr* pa[2] = { expr_new_integer(n), expr_new_integer(-1) };
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
    Expr* ta[2] = { sum, inv };
    return evaluate(expr_new_function(expr_new_symbol(SYM_Times), ta, 2));
}
