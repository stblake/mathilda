/* edgeconnectivity.c - EdgeConnectivity[g]: the minimum number of edges whose
 * removal disconnects g (0 if g is already disconnected / not strongly
 * connected). The edge companion to VertexConnectivity.
 *
 * By the max-flow/min-cut theorem the minimum s-t edge cut equals the max s-t
 * flow with unit edge capacities. Global edge connectivity is then:
 *   undirected: fix a source s0 and take min over all t != s0 of maxflow(s0,t)
 *               (a standard result — one source suffices),
 *   directed:   min over all ordered pairs (s,t), s != t.
 * Each flow is Edmonds-Karp (BFS-augmenting), O(V*E^2); overall O(V^2) flows for
 * the directed case, O(V) for undirected — fine for the MVP's small graphs.
 *
 * Memory (SPEC section 4): returns a fresh integer; frees res. NULL on a
 * non-graph argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>

/* Max s->t flow over the unit-capacity residual graph `cap` (n x n). */
static int max_flow(const int* cap, int n, int s, int t) {
    int* res = malloc((size_t)n * n * sizeof(int));
    int* parent = malloc((size_t)n * sizeof(int));
    int* q = malloc((size_t)n * sizeof(int));
    if (!res || !parent || !q) { free(res); free(parent); free(q); return 0; }
    memcpy(res, cap, (size_t)n * n * sizeof(int));

    int flow = 0;
    for (;;) {
        for (int i = 0; i < n; i++) parent[i] = -1;
        parent[s] = s;
        int head = 0, tail = 0;
        q[tail++] = s;
        while (head < tail && parent[t] < 0) {
            int u = q[head++];
            for (int v = 0; v < n; v++)
                if (parent[v] < 0 && res[(size_t)u * n + v] > 0) { parent[v] = u; q[tail++] = v; }
        }
        if (parent[t] < 0) break;                 /* no augmenting path */
        int bottleneck = INT_MAX;
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            if (res[(size_t)u * n + v] < bottleneck) bottleneck = res[(size_t)u * n + v];
        }
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            res[(size_t)u * n + v] -= bottleneck;
            res[(size_t)v * n + u] += bottleneck;
        }
        flow += bottleneck;
    }
    free(res); free(parent); free(q);
    return flow;
}

Expr* builtin_edge_connectivity(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;
    if (n < 2) return expr_new_integer(0);

    int directed = (ne > 0);
    for (size_t i = 0; i < ne; i++)
        if (graph_edge_kind(edges->data.function.args[i]) != SYM_DirectedEdge) { directed = 0; break; }

    int* cap = calloc((size_t)n * n, sizeof(int));
    if (!cap) return NULL;
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        int a = graph_vertex_index(verts, e->data.function.args[0]);
        int b = graph_vertex_index(verts, e->data.function.args[1]);
        if (a < 0 || b < 0) continue;
        cap[(size_t)a * n + b] += 1;
        if (!directed) cap[(size_t)b * n + a] += 1;
    }

    int lambda = INT_MAX;
    if (directed) {
        for (int s = 0; s < n && lambda > 0; s++)
            for (int t = 0; t < n && lambda > 0; t++)
                if (s != t) { int f = max_flow(cap, n, s, t); if (f < lambda) lambda = f; }
    } else {
        for (int t = 1; t < n && lambda > 0; t++) {
            int f = max_flow(cap, n, 0, t);
            if (f < lambda) lambda = f;
        }
    }
    free(cap);
    return expr_new_integer(lambda == INT_MAX ? 0 : lambda);
}
