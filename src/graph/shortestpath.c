/* shortestpath.c - FindShortestPath[g,s,t] and GraphDistance[g,s,t].
 *
 * Unweighted breadth-first search over the successor adjacency (out[]): for a
 * directed graph this follows edge direction; for an undirected graph out[] is
 * symmetric, so it is an ordinary shortest path. Wolfram's naming split is
 * kept: FindShortestPath returns the vertex path, GraphDistance the length.
 *
 * Unreachable target: FindShortestPath -> {} (empty list), GraphDistance ->
 * Infinity.
 *
 * Memory (SPEC section 4): returns freshly-allocated results; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* BFS from src over out[]; fills parent[] (-1 = root/unvisited) and dist[]
 * (-1 = unreached). Caller allocates parent/dist of length a->n. */
static void bfs(const GraphAdj* a, int src, int* parent, int* dist) {
    for (int i = 0; i < a->n; i++) { parent[i] = -1; dist[i] = -1; }
    int* q = calloc((size_t)(a->n > 0 ? a->n : 1), sizeof(int));
    int head = 0, tail = 0;
    dist[src] = 0; q[tail++] = src;
    while (head < tail) {
        int u = q[head++];
        for (int j = 0; j < a->outdeg[u]; j++) {
            int w = a->out[u][j];
            if (dist[w] < 0) { dist[w] = dist[u] + 1; parent[w] = u; q[tail++] = w; }
        }
    }
    free(q);
}

/* Resolve g, s, t to a GraphAdj and endpoint indices. Returns adj (caller frees)
 * or NULL; on success sets *is,*it. */
static GraphAdj* resolve(Expr* res, int* is, int* it) {
    if (res->data.function.arg_count != 3) return NULL;
    const Expr* g = res->data.function.args[0];
    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    *is = graph_vertex_index(a->verts, res->data.function.args[1]);
    *it = graph_vertex_index(a->verts, res->data.function.args[2]);
    if (*is < 0 || *it < 0) { graph_adj_free(a); return NULL; }
    return a;
}

Expr* builtin_find_shortest_path(Expr* res) {
    int is, it;
    GraphAdj* a = resolve(res, &is, &it);
    if (!a) return NULL;

    int* parent = calloc((size_t)a->n, sizeof(int));
    int* dist   = calloc((size_t)a->n, sizeof(int));
    bfs(a, is, parent, dist);

    Expr* out;
    if (dist[it] < 0) {
        out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);  /* {} */
    } else {
        int len = dist[it] + 1;
        Expr** path = calloc((size_t)len, sizeof(Expr*));
        int v = it;
        for (int k = len - 1; k >= 0; k--) { path[k] = expr_copy(a->verts->data.function.args[v]); v = parent[v]; }
        out = expr_new_function(expr_new_symbol(SYM_List), path, (size_t)len);
        free(path);
    }
    free(parent); free(dist); graph_adj_free(a);
    return out;
}

Expr* builtin_graph_distance(Expr* res) {
    int is, it;
    GraphAdj* a = resolve(res, &is, &it);
    if (!a) return NULL;

    int* parent = calloc((size_t)a->n, sizeof(int));
    int* dist   = calloc((size_t)a->n, sizeof(int));
    bfs(a, is, parent, dist);

    Expr* out = (dist[it] < 0) ? expr_new_symbol(SYM_Infinity)
                               : expr_new_integer(dist[it]);
    free(parent); free(dist); graph_adj_free(a);
    return out;
}
