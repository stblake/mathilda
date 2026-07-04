/* stronglyconnectedq.c - StronglyConnectedGraphQ[g]: True iff g is strongly
 * connected -- every vertex is reachable from every other following edge
 * directions. (For an undirected graph this coincides with ordinary
 * connectivity, since out- and in-neighbourhoods are the same.)
 *
 * A graph is strongly connected iff, from any single vertex, every other vertex
 * is reachable forward (following out-edges) AND that vertex is reachable from
 * every other, i.e. every vertex reaches it -- which is a forward reach on the
 * reverse graph (following in-edges). So two BFS from vertex 0, one over out[]
 * and one over in[], each reaching all n vertices, settle it in O(V+E). A single
 * vertex is trivially strongly connected; the empty graph is not.
 *
 * Memory (SPEC section 4): returns True/False; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Count vertices reachable from src using successor lists nbr/deg. */
static int reach_count(int n, int src, int* const* nbr, const int* deg, int* q, char* seen) {
    for (int i = 0; i < n; i++) seen[i] = 0;
    int head = 0, tail = 0, cnt = 0;
    seen[src] = 1; q[tail++] = src;
    while (head < tail) {
        int v = q[head++]; cnt++;
        for (int j = 0; j < deg[v]; j++) {
            int u = nbr[v][j];
            if (!seen[u]) { seen[u] = 1; q[tail++] = u; }
        }
    }
    return cnt;
}

Expr* builtin_strongly_connected_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GraphAdj* a = graph_build_adj(res->data.function.args[0]);
    if (!a) return NULL;
    int n = a->n;
    if (n == 0) { graph_adj_free(a); return expr_new_symbol(SYM_False); }

    int*  q    = malloc((size_t)n * sizeof(int));
    char* seen = malloc((size_t)n);
    int fwd = reach_count(n, 0, a->out, a->outdeg, q, seen);  /* 0 reaches all? */
    int bwd = reach_count(n, 0, a->in,  a->indeg,  q, seen);  /* all reach 0?   */
    free(q); free(seen);
    graph_adj_free(a);

    int strong = (fwd == n && bwd == n);
    return expr_new_symbol(strong ? SYM_True : SYM_False);
}
