/* vertexcoreness.c - VertexCoreness[g]: for each vertex, its coreness (core
 * number) -- the largest k such that the vertex still belongs to the k-core (the
 * maximal subgraph of minimum degree k).
 *
 * Batagelj-Zaversnik peeling: repeatedly remove a vertex of minimum current
 * degree; the running maximum of those removal-degrees is the core level, and
 * each removed vertex is assigned the current level. Edge direction is ignored
 * (coreness is defined on the underlying undirected graph). Results are exact
 * integers in vertex order. O(V^2) with a linear min-degree scan (small-graph
 * scale). Complements KCoreComponents (a vertex has coreness >= k iff it appears
 * in KCoreComponents[g, k]).
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_vertex_coreness(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    char* adj = (n > 0) ? calloc((size_t)n * (size_t)n, 1) : NULL;
    int*  deg = (n > 0) ? calloc((size_t)n, sizeof(int)) : NULL;
    for (size_t e = 0; e < ne; e++) {
        const Expr* ed = edges->data.function.args[e];
        int ia = graph_vertex_index(verts, ed->data.function.args[0]);
        int ib = graph_vertex_index(verts, ed->data.function.args[1]);
        if (ia < 0 || ib < 0 || ia == ib) continue;
        if (!adj[(size_t)ia * n + ib]) {
            adj[(size_t)ia * n + ib] = adj[(size_t)ib * n + ia] = 1;
            deg[ia]++; deg[ib]++;
        }
    }

    int*  core    = (n > 0) ? calloc((size_t)n, sizeof(int)) : NULL;
    char* removed = (n > 0) ? calloc((size_t)n, 1) : NULL;
    int level = 0;
    for (int step = 0; step < n; step++) {
        int v = -1, best = 0;
        for (int i = 0; i < n; i++)
            if (!removed[i] && (v < 0 || deg[i] < best)) { v = i; best = deg[i]; }
        if (best > level) level = best;         /* running max = core level */
        core[v] = level;
        removed[v] = 1;
        for (int u = 0; u < n; u++)
            if (!removed[u] && adj[(size_t)v * n + u]) deg[u]--;
    }
    free(adj); free(deg); free(removed);

    Expr** items = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    for (int i = 0; i < n; i++) items[i] = expr_new_integer(core[i]);
    free(core);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
    free(items);
    return out;
}
