/* acyclic.c - AcyclicGraphQ[g] and TopologicalSort[g].
 *
 * Both rest on Kahn's algorithm (topological sort by repeatedly removing
 * in-degree-0 vertices): a directed graph is acyclic iff every vertex is
 * eventually removed. O(V+E) over the shared integer adjacency.
 *
 *   TopologicalSort[g]  -> a vertex order with every edge pointing forward, or
 *                          $Failed if g has a directed cycle (or undirected
 *                          edges, which act as 2-cycles). Directed graphs only,
 *                          matching the Wolfram Language.
 *   AcyclicGraphQ[g]    -> True iff g has no cycle. For an all-directed graph
 *                          that means a DAG (Kahn removes all vertices); for a
 *                          graph with undirected edges it means the underlying
 *                          undirected graph is a forest (|E| = |V| - #components).
 *
 * Memory (SPEC section 4): returns freshly-allocated results; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Kahn's algorithm. Returns a malloc'd array of vertex indices in topological
 * order (caller frees) and writes the number ordered to *count. When *count
 * equals a->n the graph is a DAG; otherwise a cycle blocked the remainder. */
static int* topo_order(const GraphAdj* a, int* count) {
    int n = a->n;
    int* indeg = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* order = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* q     = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!indeg || !order || !q) { free(indeg); free(order); free(q); *count = -1; return NULL; }
    for (int i = 0; i < n; i++) indeg[i] = a->indeg[i];
    int head = 0, tail = 0;
    for (int i = 0; i < n; i++) if (indeg[i] == 0) q[tail++] = i;
    int k = 0;
    while (head < tail) {
        int u = q[head++];
        order[k++] = u;
        for (int j = 0; j < a->outdeg[u]; j++) {
            int w = a->out[u][j];
            if (--indeg[w] == 0) q[tail++] = w;
        }
    }
    free(indeg); free(q);
    *count = k;
    return order;
}

Expr* builtin_topological_sort(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GraphAdj* a = graph_build_adj(res->data.function.args[0]);
    if (!a) return NULL;
    int cnt = 0;
    int* order = topo_order(a, &cnt);
    if (!order) { graph_adj_free(a); return NULL; }

    Expr* out;
    if (cnt == a->n) {
        Expr** items = (a->n > 0) ? calloc((size_t)a->n, sizeof(Expr*)) : NULL;
        for (int i = 0; i < a->n; i++)
            items[i] = expr_copy(a->verts->data.function.args[order[i]]);
        out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)a->n);
        free(items);
    } else {
        out = expr_new_symbol(SYM_DollarFailed);   /* cyclic / undirected */
    }
    free(order); graph_adj_free(a);
    return out;
}

Expr* builtin_acyclic_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return expr_new_symbol(SYM_False);

    const Expr* edges = g->data.function.args[1];
    size_t ne = edges->data.function.arg_count;
    int all_directed = 1;
    for (size_t i = 0; i < ne; i++)
        if (graph_edge_kind(edges->data.function.args[i]) != SYM_DirectedEdge) { all_directed = 0; break; }

    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;

    int acyclic;
    if (all_directed) {
        int cnt = 0;
        int* order = topo_order(a, &cnt);
        if (!order) { graph_adj_free(a); return NULL; }
        acyclic = (cnt == a->n);          /* DAG iff every vertex removed */
        free(order);
    } else {
        /* Underlying undirected graph is a forest iff E = V - components. */
        int comps = graph_count_components(a, NULL, NULL);
        acyclic = ((int)ne == a->n - comps);
    }
    graph_adj_free(a);
    return expr_new_symbol(acyclic ? SYM_True : SYM_False);
}
