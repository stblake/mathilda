/* completegraphq.c - CompleteGraphQ[g]: True iff g is complete -- every pair of
 * distinct vertices is adjacent (the underlying undirected graph is K_n).
 *
 * Build a boolean adjacency (direction ignored, a->b or b->a both count) and
 * check that all n(n-1)/2 unordered pairs are present. Equivalently, the number
 * of distinct undirected edges equals n(n-1)/2. O(V^2). A graph with <= 1 vertex
 * is vacuously complete (K_0, K_1).
 *
 * Memory (SPEC section 4): returns True/False; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_complete_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    char* adj = (n > 0) ? calloc((size_t)n * (size_t)n, 1) : NULL;
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        int ia = graph_vertex_index(verts, e->data.function.args[0]);
        int ib = graph_vertex_index(verts, e->data.function.args[1]);
        if (ia < 0 || ib < 0 || ia == ib) continue;
        adj[(size_t)ia * n + ib] = adj[(size_t)ib * n + ia] = 1;
    }

    int complete = 1;
    for (int i = 0; i < n && complete; i++)
        for (int j = i + 1; j < n; j++)
            if (!adj[(size_t)i * n + j]) { complete = 0; break; }
    free(adj);
    return expr_new_symbol(complete ? SYM_True : SYM_False);
}
