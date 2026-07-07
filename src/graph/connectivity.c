/* connectivity.c - ConnectedGraphQ[g] and VertexConnectivity[g].
 *
 * Both operate on the underlying undirected graph.
 *   ConnectedGraphQ[g]    True iff g has >= 1 vertex and forms a single
 *                         connected component.
 *   VertexConnectivity[g] the least number of vertices whose removal
 *                         disconnects g (n-1 for a complete graph, 0 if already
 *                         disconnected or trivial). Computed by brute-force
 *                         search over vertex subsets -- exact, intended for the
 *                         small graphs of a pico-CAS.
 *
 * Memory (SPEC section 4): returns freshly-allocated results; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_connected_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GraphAdj* a = graph_build_adj(res->data.function.args[0]);
    if (!a) return NULL;
    int comps = graph_count_components(a, NULL, NULL);
    int connected = (a->n >= 1 && comps == 1);
    graph_adj_free(a);
    return expr_new_symbol(connected ? SYM_True : SYM_False);
}

/* Does some size-k subset of vertices disconnect the graph? Enumerates all
 * C(n,k) subsets via an index array; sets `removed` for each. */
static int some_cut_of_size(const GraphAdj* a, int k, char* removed) {
    int n = a->n;
    int* idx = calloc((size_t)(k > 0 ? k : 1), sizeof(int));
    for (int i = 0; i < k; i++) idx[i] = i;
    int found = 0;
    while (!found) {
        for (int i = 0; i < n; i++) removed[i] = 0;
        for (int i = 0; i < k; i++) removed[idx[i]] = 1;
        int active = 0;
        int comps = graph_count_components(a, removed, &active);
        if (active >= 2 && comps >= 2) { found = 1; break; }
        /* Advance the combination (rightmost index that can move). */
        int p = k - 1;
        while (p >= 0 && idx[p] == n - k + p) p--;
        if (p < 0) break;
        idx[p]++;
        for (int i = p + 1; i < k; i++) idx[i] = idx[i - 1] + 1;
    }
    free(idx);
    return found;
}

Expr* builtin_vertex_connectivity(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GraphAdj* a = graph_build_adj(res->data.function.args[0]);
    if (!a) return NULL;
    int n = a->n;

    int kappa;
    if (n <= 1) {
        kappa = 0;
    } else if (graph_count_components(a, NULL, NULL) > 1) {
        kappa = 0;                       /* already disconnected */
    } else {
        char* removed = calloc((size_t)n, sizeof(char));
        kappa = n - 1;                   /* complete-graph fallback */
        for (int k = 1; k <= n - 2; k++) {
            if (some_cut_of_size(a, k, removed)) { kappa = k; break; }
        }
        free(removed);
    }
    graph_adj_free(a);
    return expr_new_integer(kappa);
}
