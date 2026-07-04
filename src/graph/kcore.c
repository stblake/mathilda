/* kcore.c - KCoreComponents[g, k]: the connected components of the k-core of g,
 * as a list of vertex lists. The k-core is the maximal subgraph in which every
 * vertex has (undirected) degree at least k; it is found by repeatedly deleting
 * any vertex whose current degree drops below k until none remain.
 *
 * Peeling is O(V + E): each vertex is queued for removal at most once, and each
 * edge is examined a constant number of times as its endpoints are deleted. The
 * survivors are then split into connected components by BFS. Edge direction is
 * ignored (the k-core is defined on the underlying undirected graph), matching
 * Wolfram; components come out ordered by least vertex index, vertices within a
 * component in canonical order.
 *
 * k must be a non-negative integer literal, else the call is left unevaluated.
 *
 * Memory (SPEC section 4): returns a fresh List of Lists; frees res. NULL on a
 * non-graph argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_kcore_components(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g  = res->data.function.args[0];
    const Expr* ke = res->data.function.args[1];
    if (!graph_is_valid(g)) return NULL;
    if (ke->type != EXPR_INTEGER || ke->data.integer < 0) return NULL;
    long k = (long)ke->data.integer;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    /* Symmetric boolean adjacency (direction ignored) + degrees. */
    char* adj = (n > 0) ? calloc((size_t)n * (size_t)n, 1) : NULL;
    int*  deg = (n > 0) ? calloc((size_t)n, sizeof(int)) : NULL;
    for (size_t e = 0; e < ne; e++) {
        const Expr* ed = edges->data.function.args[e];
        int ia = graph_vertex_index(verts, ed->data.function.args[0]);
        int ib = graph_vertex_index(verts, ed->data.function.args[1]);
        if (ia < 0 || ib < 0 || ia == ib) continue;
        if (!adj[(size_t)ia * n + ib]) {           /* first time this pair */
            adj[(size_t)ia * n + ib] = adj[(size_t)ib * n + ia] = 1;
            deg[ia]++; deg[ib]++;
        }
    }

    /* Peel: remove any alive vertex whose degree < k, cascading. */
    char* alive = (n > 0) ? malloc((size_t)n) : NULL;
    for (int i = 0; i < n; i++) alive[i] = 1;
    int* queue = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;
    int qh = 0, qt = 0;
    for (int i = 0; i < n; i++) if (deg[i] < k) queue[qt++] = i;
    while (qh < qt) {
        int v = queue[qh++];
        if (!alive[v]) continue;
        alive[v] = 0;
        for (int u = 0; u < n; u++)
            if (adj[(size_t)v * n + u] && alive[u] && --deg[u] == k - 1)
                queue[qt++] = u;   /* enqueue once, on the downward crossing */
    }

    /* BFS the survivors into components, ordered by least vertex index. */
    char* seen = (n > 0) ? calloc((size_t)n, 1) : NULL;
    Expr** comps = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    size_t ncomp = 0;
    for (int s = 0; s < n; s++) {
        if (!alive[s] || seen[s]) continue;
        int bh = 0, bt = 0;
        queue[bt++] = s; seen[s] = 1;
        int* members = malloc((size_t)n * sizeof(int));
        int mc = 0;
        while (bh < bt) {
            int v = queue[bh++];
            members[mc++] = v;
            for (int u = 0; u < n; u++)
                if (adj[(size_t)v * n + u] && alive[u] && !seen[u]) { seen[u] = 1; queue[bt++] = u; }
        }
        /* members collected in BFS order; emit in canonical (index) order. */
        Expr** vs = calloc((size_t)mc, sizeof(Expr*));
        int idx = 0;
        for (int v = 0; v < n; v++) {
            int in_comp = 0;
            for (int m = 0; m < mc; m++) if (members[m] == v) { in_comp = 1; break; }
            if (in_comp) vs[idx++] = expr_copy(verts->data.function.args[v]);
        }
        comps[ncomp++] = expr_new_function(expr_new_symbol(SYM_List), vs, (size_t)mc);
        free(vs); free(members);
    }

    free(adj); free(deg); free(alive); free(queue); free(seen);

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), comps, ncomp);
    free(comps);
    return out;
}
