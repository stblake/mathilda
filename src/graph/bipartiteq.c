/* bipartiteq.c - BipartiteGraphQ[g]: True iff the underlying undirected graph
 * is 2-colorable (no odd cycle), i.e. its vertices split into two independent
 * sets. Edge direction is ignored (a graph is bipartite as an undirected
 * structure). Vacuously True for an edgeless graph.
 *
 * Algorithm: a single BFS per component, 2-coloring as it goes; a bipartite
 * graph is exactly one with no edge joining two same-colored vertices. O(V+E)
 * over the on-demand integer adjacency (graph_build_adj), which stores each
 * undirected edge symmetrically in out[]/in[], so scanning both neighbor lists
 * visits the whole undirected neighborhood.
 *
 * Memory (SPEC section 4): returns a fresh True/False symbol; the evaluator
 * frees res. False (not NULL) on a non-graph argument, matching GraphQ.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_bipartite_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return expr_new_symbol(SYM_False);

    GraphAdj* a = graph_build_adj(g);
    if (!a) return expr_new_symbol(SYM_False);
    int n = a->n;

    int bipartite = 1;
    int* color = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;
    int* queue = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;
    if (n > 0 && (!color || !queue)) { free(color); free(queue); graph_adj_free(a); return NULL; }
    for (int i = 0; i < n; i++) color[i] = -1;

    for (int s = 0; s < n && bipartite; s++) {
        if (color[s] >= 0) continue;
        int head = 0, tail = 0;
        color[s] = 0; queue[tail++] = s;
        while (head < tail && bipartite) {
            int u = queue[head++];
            /* Visit both stored directions to cover the undirected neighborhood. */
            for (int pass = 0; pass < 2 && bipartite; pass++) {
                int deg = pass ? a->indeg[u] : a->outdeg[u];
                int* nb = pass ? a->in[u] : a->out[u];
                for (int e = 0; e < deg; e++) {
                    int v = nb[e];
                    if (color[v] < 0) { color[v] = 1 - color[u]; queue[tail++] = v; }
                    else if (color[v] == color[u]) { bipartite = 0; break; }
                }
            }
        }
    }

    free(color); free(queue);
    graph_adj_free(a);
    return expr_new_symbol(bipartite ? SYM_True : SYM_False);
}
