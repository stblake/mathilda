/* neighborhoodgraph.c - NeighborhoodGraph[g, v] / NeighborhoodGraph[g, v, k]: the
 * subgraph of g induced by v together with every vertex within graph distance k
 * of v (k defaults to 1). The induced subgraph keeps all edges of g whose both
 * endpoints lie in that neighbourhood.
 *
 * A depth-limited BFS from v over the successor adjacency (out[]) collects the
 * neighbourhood -- direction-aware, matching GraphDistance -- then the edges of g
 * between kept vertices are retained (kinds preserved). O(V + E). k must be a
 * non-negative integer; k = 0 gives the single vertex v.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL on a
 * non-graph, an unknown vertex, or a negative/non-integer k.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_neighborhood_graph(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 2 && argc != 3) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    long k = 1;
    if (argc == 3) {
        const Expr* ke = res->data.function.args[2];
        if (ke->type != EXPR_INTEGER || ke->data.integer < 0) return NULL;
        k = (long)ke->data.integer;
    }

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    int src = graph_vertex_index(verts, res->data.function.args[1]);
    if (src < 0) return NULL;                        /* v is not a vertex */

    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;

    /* Depth-limited BFS from src: keep[i] set for dist(src, i) <= k. */
    char* keep = calloc((size_t)n, 1);
    int* dist  = malloc((size_t)n * sizeof(int));
    int* q     = malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) dist[i] = -1;
    int head = 0, tail = 0;
    dist[src] = 0; keep[src] = 1; q[tail++] = src;
    while (head < tail) {
        int v = q[head++];
        if (dist[v] >= (int)k) continue;
        for (int j = 0; j < a->outdeg[v]; j++) {
            int u = a->out[v][j];
            if (dist[u] < 0) { dist[u] = dist[v] + 1; keep[u] = 1; q[tail++] = u; }
        }
    }
    free(dist); free(q);
    graph_adj_free(a);

    /* Induced subgraph on the kept vertices. */
    Expr** vc = malloc((size_t)(n > 0 ? n : 1) * sizeof(Expr*));
    size_t nv = 0;
    for (int i = 0; i < n; i++) if (keep[i]) vc[nv++] = expr_copy(verts->data.function.args[i]);

    size_t m = edges->data.function.arg_count;
    Expr** ec = (m > 0) ? malloc(m * sizeof(Expr*)) : NULL;
    size_t me = 0;
    for (size_t e = 0; e < m; e++) {
        const Expr* ed = edges->data.function.args[e];
        int ia = graph_vertex_index(verts, ed->data.function.args[0]);
        int ib = graph_vertex_index(verts, ed->data.function.args[1]);
        if (ia >= 0 && ib >= 0 && keep[ia] && keep[ib]) ec[me++] = expr_copy((Expr*)ed);
    }
    free(keep);

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, nv);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
