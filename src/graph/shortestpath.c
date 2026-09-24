/* shortestpath.c - FindShortestPath[g,s,t] and GraphDistance[g,s,t].
 *
 * Two algorithms, dispatched on graph_weights_usable(g):
 *   - Unweighted (default): breadth-first search over the successor adjacency
 *     (GraphAdj.out[]): for a directed graph this follows edge direction; for
 *     an undirected graph out[] is symmetric, so it is an ordinary shortest
 *     path.
 *   - Weighted (g carries a non-negative-numeric EdgeWeight): Dijkstra over a
 *     local, call-scoped weighted adjacency (WAdj below) -- NOT over
 *     GraphAdj, which has no weight storage and is shared by 5 other
 *     builtins (ConnectedComponents, WeaklyConnectedComponents,
 *     FindSpanningTree, ConnectedGraphQ, VertexConnectivity); widening it
 *     would risk the exact class of shared-choke-point defect a
 *     plan-reviewer pass caught during the EdgeWeight ticket. Falls back to
 *     BFS for a symbolic or negative weight rather than erroring.
 *
 * Wolfram's naming split is kept: FindShortestPath returns the vertex path,
 * GraphDistance the length/total weight.
 *
 * Unreachable target: FindShortestPath -> {} (empty list), GraphDistance ->
 * Infinity.
 *
 * Memory (SPEC section 4): returns freshly-allocated results; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include <stdlib.h>
#include <float.h>

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

/* ---- Weighted adjacency, local to this file (see header comment) --------- */

typedef struct {
    int n;
    const Expr* verts;   /* borrowed from g */
    int*    outdeg;
    int**   out;         /* out[v][k]   = neighbor vertex index */
    Expr*** outw;        /* outw[v][k]  = borrowed weight Expr* for that edge */
    Expr*   weights;     /* owns every Expr* referenced by outw; freed by wadj_free */
} WAdj;

static void wadj_free(WAdj* w) {
    if (!w) return;
    for (int i = 0; i < w->n; i++) { free(w->out[i]); free(w->outw[i]); }
    free(w->out); free(w->outw); free(w->outdeg);
    if (w->weights) expr_free(w->weights);
    free(w);
}

/* Builds a weighted successor adjacency directly from g's canonical form
 * (verts/edges/EdgeWeight), independent of GraphAdj. NULL if g is not a
 * valid, weights-usable graph. */
static WAdj* build_wadj(const Expr* g) {
    if (!graph_weights_usable(g)) return NULL;
    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    Expr* weights = graph_resolve_edge_weights(g);
    if (!weights) return NULL;

    WAdj* w = calloc(1, sizeof(WAdj));
    if (!w) { expr_free(weights); return NULL; }
    w->n = n; w->verts = verts; w->weights = weights;
    w->outdeg = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    w->out    = calloc((size_t)(n > 0 ? n : 1), sizeof(int*));
    w->outw   = calloc((size_t)(n > 0 ? n : 1), sizeof(Expr**));

    GraphVIdx* ix = graph_vidx_new((size_t)n);
    for (int i = 0; i < n; i++)
        graph_vidx_put(ix, verts->data.function.args[i], i);

    /* Pass 1: count out-degrees (each undirected edge contributes to both endpoints). */
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        const char* kind = graph_edge_kind(e);
        int ia = graph_vidx_get(ix, e->data.function.args[0]);
        int ib = graph_vidx_get(ix, e->data.function.args[1]);
        if (ia < 0 || ib < 0) continue;
        w->outdeg[ia]++;
        if (kind == SYM_UndirectedEdge) w->outdeg[ib]++;
    }
    for (int i = 0; i < n; i++) {
        w->out[i]  = w->outdeg[i] ? calloc((size_t)w->outdeg[i], sizeof(int)) : NULL;
        w->outw[i] = w->outdeg[i] ? calloc((size_t)w->outdeg[i], sizeof(Expr*)) : NULL;
    }
    int* fill = calloc((size_t)(n > 0 ? n : 1), sizeof(int));

    /* Pass 2: fill. */
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        const char* kind = graph_edge_kind(e);
        int ia = graph_vidx_get(ix, e->data.function.args[0]);
        int ib = graph_vidx_get(ix, e->data.function.args[1]);
        if (ia < 0 || ib < 0) continue;
        Expr* wt = weights->data.function.args[k];
        w->out[ia][fill[ia]]  = ib;
        w->outw[ia][fill[ia]] = wt;
        fill[ia]++;
        if (kind == SYM_UndirectedEdge) {
            w->out[ib][fill[ib]]  = ia;
            w->outw[ib][fill[ib]] = wt;
            fill[ib]++;
        }
    }
    free(fill);
    graph_vidx_free(ix);
    return w;
}

/* Dijkstra from src over w; fills parent[] (-1 = root/unvisited) and dist[]
 * (DBL_MAX = unreached). O(V^2) array scan, no heap -- consistent with this
 * subsystem's existing small-graph exact-algorithm precedent
 * (VertexConnectivity's own brute-force). dist[t] is also GraphDistance's
 * answer: a machine real, as in the Wolfram Language. */
static void dijkstra(const WAdj* w, int src, int* parent, double* dist) {
    char* done = calloc((size_t)(w->n > 0 ? w->n : 1), sizeof(char));
    for (int i = 0; i < w->n; i++) { parent[i] = -1; dist[i] = DBL_MAX; }
    dist[src] = 0.0;

    for (int iter = 0; iter < w->n; iter++) {
        int u = -1;
        double best = DBL_MAX;
        for (int i = 0; i < w->n; i++)
            if (!done[i] && dist[i] < best) { best = dist[i]; u = i; }
        if (u < 0) break;                 /* remaining vertices are unreachable */
        done[u] = 1;
        for (int j = 0; j < w->outdeg[u]; j++) {
            int v = w->out[u][j];
            double d = dist[u] + graph_weight_to_double(w->outw[u][j]);
            if (d < dist[v]) { dist[v] = d; parent[v] = u; }
        }
    }
    free(done);
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
    const Expr* g = res->data.function.args[0];

    Expr* out;
    if (graph_weights_usable(g)) {
        WAdj* w = build_wadj(g);
        int* parent = calloc((size_t)w->n, sizeof(int));
        double* dist = calloc((size_t)w->n, sizeof(double));
        dijkstra(w, is, parent, dist);
        if (dist[it] == DBL_MAX) {
            out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        } else {
            /* Reconstruct the vertex path from parent[] (same shape as BFS's). */
            int len = 1, v = it;
            while (v != is) { len++; v = parent[v]; }
            Expr** path = calloc((size_t)len, sizeof(Expr*));
            v = it;
            for (int k = len - 1; k >= 0; k--) { path[k] = expr_copy(w->verts->data.function.args[v]); v = parent[v]; }
            out = expr_new_function(expr_new_symbol(SYM_List), path, (size_t)len);
            free(path);
        }
        free(parent); free(dist); wadj_free(w);
    } else {
        int* parent = calloc((size_t)a->n, sizeof(int));
        int* dist   = calloc((size_t)a->n, sizeof(int));
        bfs(a, is, parent, dist);
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
        free(parent); free(dist);
    }
    graph_adj_free(a);
    return out;
}

Expr* builtin_graph_distance(Expr* res) {
    int is, it;
    GraphAdj* a = resolve(res, &is, &it);
    if (!a) return NULL;
    const Expr* g = res->data.function.args[0];

    Expr* out;
    if (graph_weights_usable(g)) {
        WAdj* w = build_wadj(g);
        int* parent = calloc((size_t)w->n, sizeof(int));
        double* dist = calloc((size_t)w->n, sizeof(double));
        dijkstra(w, is, parent, dist);
        /* A weighted distance is a machine real, as in the Wolfram Language
         * (GraphDistance[g, 1, 3] with weights {5, 7} is 12., not 12): the
         * Dijkstra accumulator itself, summed along the path in path order, so
         * it agrees bit-for-bit with GraphDistance[g, s] / GraphDistanceMatrix. */
        out = (dist[it] == DBL_MAX) ? expr_new_symbol(SYM_Infinity)
                                    : expr_new_real(dist[it]);
        free(parent); free(dist); wadj_free(w);
    } else {
        int* parent = calloc((size_t)a->n, sizeof(int));
        int* dist   = calloc((size_t)a->n, sizeof(int));
        bfs(a, is, parent, dist);
        out = (dist[it] < 0) ? expr_new_symbol(SYM_Infinity)
                             : expr_new_integer(dist[it]);
        free(parent); free(dist);
    }
    graph_adj_free(a);
    return out;
}
