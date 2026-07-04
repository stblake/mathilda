/* centrality.c - ClosenessCentrality[g].
 *
 * For each vertex i let r_i be the number of vertices reachable from i (itself
 * included) and S_i the sum of shortest-path distances from i to those reachable
 * vertices. Closeness centrality (Wolfram's disconnected-safe form) is
 *     c_i = (r_i - 1)^2 / ((n - 1) * S_i),
 * which for a connected graph reduces to the usual (n - 1) / S_i, and is 0 for an
 * isolated vertex. Distances follow edge direction on a directed graph.
 *
 * The result is exact: each c_i is built as num/den and reduced by the evaluator,
 * so a symmetric graph yields clean rationals (e.g. every cycle vertex is 3/4 on
 * C4). O(V*(V+E)) — a BFS per vertex over the shared adjacency.
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include "eval.h"
#include <stdlib.h>

/* BFS from src over out[]; dist[i] = distance or -1 if unreachable. */
static void bfs_out(const GraphAdj* a, int src, int* dist, int* q) {
    int n = a->n;
    for (int i = 0; i < n; i++) dist[i] = -1;
    int head = 0, tail = 0;
    dist[src] = 0; q[tail++] = src;
    while (head < tail) {
        int u = q[head++];
        for (int j = 0; j < a->outdeg[u]; j++) {
            int w = a->out[u][j];
            if (dist[w] < 0) { dist[w] = dist[u] + 1; q[tail++] = w; }
        }
    }
}

/* Exact num/den as a reduced value (integer or Rational). den > 0 required. */
static Expr* ratio(long num, long den) {
    Expr* pa[2] = { expr_new_integer(den), expr_new_integer(-1) };
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
    Expr* ta[2] = { expr_new_integer(num), inv };
    return evaluate(expr_new_function(expr_new_symbol(SYM_Times), ta, 2));
}

Expr* builtin_closeness_centrality(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GraphAdj* a = graph_build_adj(res->data.function.args[0]);
    if (!a) return NULL;
    int n = a->n;

    int* dist = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* q    = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    Expr** items = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    if ((n > 0 && (!dist || !q || !items))) { free(dist); free(q); free(items); graph_adj_free(a); return NULL; }

    for (int i = 0; i < n; i++) {
        bfs_out(a, i, dist, q);
        long reach = 0, sum = 0;
        for (int j = 0; j < n; j++) if (dist[j] >= 0) { reach++; sum += dist[j]; }
        long den = (long)(n - 1) * sum;             /* (n-1) * S_i */
        items[i] = (den <= 0) ? expr_new_integer(0)
                              : ratio((reach - 1) * (reach - 1), den);
    }
    free(dist); free(q);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
    free(items);
    graph_adj_free(a);
    return out;
}
