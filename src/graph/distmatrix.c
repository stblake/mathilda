/* distmatrix.c - GraphDistanceMatrix[g]: the n x n matrix whose (i,j) entry is
 * the shortest-path distance from vertex i to vertex j -- 0 on the diagonal and
 * Infinity when j is unreachable from i. The all-pairs companion to
 * GraphDistance[g, s, t] and the distance metrics in metrics.c.
 *
 * One breadth-first search per source over the successor adjacency (out[]): for
 * a directed graph this follows edge direction (so the matrix is generally
 * asymmetric), for an undirected graph out[] is symmetric so it is ordinary
 * shortest-path distance. Complexity O(V*(V+E)).
 *
 * Memory (SPEC section 4): returns a freshly-built List of Lists; frees res.
 * Returns NULL (unevaluated) on a non-graph argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_graph_distance_matrix(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    int n = a->n;

    int* dist = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;
    int* q    = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;

    Expr** rows = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    for (int s = 0; s < n; s++) {
        for (int i = 0; i < n; i++) dist[i] = -1;
        int head = 0, tail = 0;
        dist[s] = 0; q[tail++] = s;
        while (head < tail) {
            int u = q[head++];
            for (int j = 0; j < a->outdeg[u]; j++) {
                int w = a->out[u][j];
                if (dist[w] < 0) { dist[w] = dist[u] + 1; q[tail++] = w; }
            }
        }
        Expr** entries = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
        for (int t = 0; t < n; t++)
            entries[t] = (dist[t] < 0) ? expr_new_symbol(SYM_Infinity)
                                       : expr_new_integer(dist[t]);
        rows[s] = expr_new_function(expr_new_symbol(SYM_List), entries, (size_t)n);
        free(entries);
    }
    free(dist); free(q);
    graph_adj_free(a);

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    return out;
}
