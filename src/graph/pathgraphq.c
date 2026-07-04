/* pathgraphq.c - PathGraphQ[g]: True iff g is a path graph -- the vertices can be
 * arranged in a line so that consecutive ones are joined and no others are.
 *
 * A path is exactly a tree whose maximum degree is at most 2, so the test is:
 * connected (on the underlying undirected graph) with n-1 distinct edges (tree)
 * and every vertex of undirected degree <= 2. Edge direction is ignored, and
 * a->b together with b->a counts as one undirected edge. O(V^2) via the boolean
 * adjacency. A single vertex is a (trivial) path; the empty graph and any
 * disconnected, cyclic, or branching graph is not.
 *
 * Memory (SPEC section 4): returns True/False; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_path_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;
    if (n == 0) return expr_new_symbol(SYM_False);

    char* adj = calloc((size_t)n * (size_t)n, 1);
    int*  deg = calloc((size_t)n, sizeof(int));
    long m = 0;
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        int ia = graph_vertex_index(verts, e->data.function.args[0]);
        int ib = graph_vertex_index(verts, e->data.function.args[1]);
        if (ia < 0 || ib < 0 || ia == ib) continue;
        if (!adj[(size_t)ia * n + ib]) {
            adj[(size_t)ia * n + ib] = adj[(size_t)ib * n + ia] = 1;
            deg[ia]++; deg[ib]++; m++;
        }
    }

    /* Max degree <= 2 (a path never branches). */
    int ok = 1;
    for (int i = 0; i < n && ok; i++) if (deg[i] > 2) ok = 0;

    /* Connected via BFS over the undirected adjacency. */
    int reached = 0;
    if (ok) {
        char* seen = calloc((size_t)n, 1);
        int* q = malloc((size_t)n * sizeof(int));
        int head = 0, tail = 0;
        seen[0] = 1; q[tail++] = 0;
        while (head < tail) {
            int v = q[head++]; reached++;
            for (int u = 0; u < n; u++)
                if (adj[(size_t)v * n + u] && !seen[u]) { seen[u] = 1; q[tail++] = u; }
        }
        free(seen); free(q);
    }
    free(adj); free(deg);

    int is_path = ok && (reached == n) && (m == (long)n - 1);   /* tree + deg<=2 */
    return expr_new_symbol(is_path ? SYM_True : SYM_False);
}
