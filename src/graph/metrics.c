/* metrics.c - distance metrics: VertexEccentricity, GraphDiameter,
 * GraphRadius, GraphCenter.
 *
 * All derive from single-source BFS over the successor adjacency (out[]): for a
 * directed graph this follows edge direction, for an undirected graph out[] is
 * symmetric so it is ordinary shortest-path distance. The eccentricity of a
 * vertex is the greatest distance from it to any other vertex, or Infinity if
 * some vertex is unreachable (possible per-vertex for directed graphs). Then:
 *   diameter = max eccentricity  (Infinity if any vertex has infinite ecc),
 *   radius   = min finite eccentricity  (Infinity if none is finite),
 *   center   = the vertices attaining the radius.
 *
 * Complexity O(V*(V+E)) — a BFS per vertex over the on-demand integer adjacency.
 *
 * Memory (SPEC section 4): returns freshly-allocated results; frees res. Returns
 * NULL (unevaluated) on a non-graph argument or an unknown vertex.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Eccentricity of vertex `src`: max BFS distance to any vertex, or -1 if some
 * vertex is unreachable. */
static int eccentricity(const GraphAdj* a, int src, int* dist, int* q) {
    int n = a->n;
    for (int i = 0; i < n; i++) dist[i] = -1;
    int head = 0, tail = 0, ecc = 0;
    dist[src] = 0; q[tail++] = src;
    while (head < tail) {
        int u = q[head++];
        if (dist[u] > ecc) ecc = dist[u];
        for (int j = 0; j < a->outdeg[u]; j++) {
            int w = a->out[u][j];
            if (dist[w] < 0) { dist[w] = dist[u] + 1; q[tail++] = w; }
        }
    }
    for (int i = 0; i < n; i++) if (dist[i] < 0) return -1;   /* unreachable */
    return ecc;
}

/* Compute eccentricities for every vertex into ecc[] (-1 = infinite). Returns
 * the built adjacency (caller frees) or NULL. Allocates/frees scratch itself. */
static GraphAdj* all_eccentricities(const Expr* g, int** ecc_out) {
    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    int n = a->n;
    int* ecc  = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* dist = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* q    = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!ecc || !dist || !q) { free(ecc); free(dist); free(q); graph_adj_free(a); return NULL; }
    for (int i = 0; i < n; i++) ecc[i] = eccentricity(a, i, dist, q);
    free(dist); free(q);
    *ecc_out = ecc;
    return a;
}

/* VertexEccentricity[g] -> list; VertexEccentricity[g, v] -> one value. */
Expr* builtin_vertex_eccentricity(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    int* ecc = NULL;
    GraphAdj* a = all_eccentricities(g, &ecc);
    if (!a) return NULL;
    int n = a->n;

    Expr* out;
    if (argc == 2) {
        int v = graph_vertex_index(a->verts, res->data.function.args[1]);
        if (v < 0) { free(ecc); graph_adj_free(a); return NULL; }
        out = (ecc[v] < 0) ? expr_new_symbol(SYM_Infinity) : expr_new_integer(ecc[v]);
    } else {
        Expr** items = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
        for (int i = 0; i < n; i++)
            items[i] = (ecc[i] < 0) ? expr_new_symbol(SYM_Infinity) : expr_new_integer(ecc[i]);
        out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
        free(items);
    }
    free(ecc); graph_adj_free(a);
    return out;
}

Expr* builtin_graph_diameter(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    int* ecc = NULL;
    GraphAdj* a = all_eccentricities(res->data.function.args[0], &ecc);
    if (!a) return NULL;
    int n = a->n, diam = 0, infinite = 0;
    for (int i = 0; i < n; i++) {
        if (ecc[i] < 0) { infinite = 1; break; }
        if (ecc[i] > diam) diam = ecc[i];
    }
    free(ecc); graph_adj_free(a);
    return infinite ? expr_new_symbol(SYM_Infinity) : expr_new_integer(diam);
}

Expr* builtin_graph_radius(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    int* ecc = NULL;
    GraphAdj* a = all_eccentricities(res->data.function.args[0], &ecc);
    if (!a) return NULL;
    int n = a->n, radius = -1;
    for (int i = 0; i < n; i++)
        if (ecc[i] >= 0 && (radius < 0 || ecc[i] < radius)) radius = ecc[i];
    free(ecc); graph_adj_free(a);
    return (radius < 0) ? expr_new_symbol(SYM_Infinity) : expr_new_integer(radius);
}

/* GraphCenter[g] -> vertices attaining the (finite) radius; {} if the radius is
 * Infinity (no vertex reaches all others). */
Expr* builtin_graph_center(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    int* ecc = NULL;
    GraphAdj* a = all_eccentricities(res->data.function.args[0], &ecc);
    if (!a) return NULL;
    int n = a->n, radius = -1;
    for (int i = 0; i < n; i++)
        if (ecc[i] >= 0 && (radius < 0 || ecc[i] < radius)) radius = ecc[i];

    Expr** items = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    size_t k = 0;
    if (radius >= 0)
        for (int i = 0; i < n; i++)
            if (ecc[i] == radius) items[k++] = expr_copy(a->verts->data.function.args[i]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, k);
    free(items); free(ecc); graph_adj_free(a);
    return out;
}

/* GraphPeriphery[g]: the vertices whose eccentricity equals the graph diameter.
 * When some vertex has infinite eccentricity (the graph is not strongly
 * connected), the diameter is Infinity and the periphery is exactly those
 * infinite-eccentricity vertices. The dual of GraphCenter. */
Expr* builtin_graph_periphery(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    int* ecc = NULL;
    GraphAdj* a = all_eccentricities(res->data.function.args[0], &ecc);
    if (!a) return NULL;
    int n = a->n, infinite = 0, diam = -1;
    for (int i = 0; i < n; i++) {
        if (ecc[i] < 0) infinite = 1;
        else if (ecc[i] > diam) diam = ecc[i];
    }
    Expr** items = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    size_t k = 0;
    for (int i = 0; i < n; i++) {
        int on_periphery = infinite ? (ecc[i] < 0) : (ecc[i] == diam && diam >= 0);
        if (on_periphery) items[k++] = expr_copy(a->verts->data.function.args[i]);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, k);
    free(items); free(ecc); graph_adj_free(a);
    return out;
}
