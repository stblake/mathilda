/* treegraphq.c - TreeGraphQ[g]: True iff g is a tree -- a connected graph with
 * no cycles. Equivalently: connected (on the underlying undirected graph) with
 * exactly n-1 distinct edges, for n >= 1 vertices.
 *
 * A connected graph with n-1 edges is necessarily acyclic, so the two cheap
 * checks -- underlying connectivity (one BFS) and a distinct-undirected-edge
 * count of n-1 -- together characterise a tree without a separate cycle scan.
 * Edge direction is ignored, and a->b together with b->a counts as one
 * undirected edge (a boolean adjacency dedups them). O(V^2) via the adjacency
 * matrix. A single vertex is a tree; the empty graph and any disconnected or
 * edgeless multi-vertex graph is not.
 *
 * Memory (SPEC section 4): returns True/False; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_tree_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    if (n == 0) return expr_new_symbol(SYM_False);   /* no vertices: not a tree */

    char* adj = calloc((size_t)n * (size_t)n, 1);
    long m = 0;                                       /* distinct undirected edges */
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        int ia = graph_vertex_index(verts, e->data.function.args[0]);
        int ib = graph_vertex_index(verts, e->data.function.args[1]);
        if (ia < 0 || ib < 0 || ia == ib) continue;
        if (!adj[(size_t)ia * n + ib]) {
            adj[(size_t)ia * n + ib] = adj[(size_t)ib * n + ia] = 1;
            m++;
        }
    }

    /* BFS from vertex 0 over the undirected adjacency; count reached vertices. */
    char* seen = calloc((size_t)n, 1);
    int* q = malloc((size_t)n * sizeof(int));
    int head = 0, tail = 0, reached = 0;
    seen[0] = 1; q[tail++] = 0;
    while (head < tail) {
        int v = q[head++]; reached++;
        for (int u = 0; u < n; u++)
            if (adj[(size_t)v * n + u] && !seen[u]) { seen[u] = 1; q[tail++] = u; }
    }
    free(adj); free(seen); free(q);

    int is_tree = (reached == n) && (m == (long)n - 1);
    return expr_new_symbol(is_tree ? SYM_True : SYM_False);
}
