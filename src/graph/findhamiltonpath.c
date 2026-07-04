/* findhamiltonpath.c - FindHamiltonianPath[g]: a Hamiltonian path (a walk
 * visiting every vertex exactly once, not necessarily returning to the start)
 * as a vertex list {v0, v1, ..., v_{n-1}}, or {} if none exists.
 *
 * Depth-first backtracking with visited-set pruning. Unlike a Hamiltonian
 * *cycle* -- where any vertex can serve as the root because the tour is closed
 * -- a path's endpoints are free, so the search is retried from every start
 * vertex until one yields a full-length path. Successors come from
 * GraphAdj.out[], which encodes edge direction (an undirected edge sits in both
 * endpoints' out[]), so directed and undirected graphs share one code path.
 *
 * O(V!) worst case, tiny in practice on gallery-sized graphs. Deterministic:
 * starts are tried in vertex order, neighbours in adjacency order.
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Extend path[0..depth-1]; returns 1 with path[0..n-1] filled on success. */
static int hp_extend(const GraphAdj* a, int* path, char* visited, int depth) {
    if (depth == a->n) return 1;
    int last = path[depth - 1];
    for (int j = 0; j < a->outdeg[last]; j++) {
        int u = a->out[last][j];
        if (!visited[u]) {
            visited[u] = 1;
            path[depth] = u;
            if (hp_extend(a, path, visited, depth + 1)) return 1;
            visited[u] = 0;
        }
    }
    return 0;
}

Expr* builtin_find_hamiltonian_path(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    int n = a->n;

    Expr* out = NULL;
    if (n == 0) {
        out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    } else {
        int*  path    = malloc((size_t)n * sizeof(int));
        char* visited = calloc((size_t)n, 1);
        int found = 0;
        for (int s = 0; s < n && !found; s++) {
            for (int i = 0; i < n; i++) visited[i] = 0;
            path[0] = s; visited[s] = 1;
            found = hp_extend(a, path, visited, 1);
        }
        if (found) {
            Expr** items = calloc((size_t)n, sizeof(Expr*));
            for (int i = 0; i < n; i++)
                items[i] = expr_copy((Expr*)a->verts->data.function.args[path[i]]);
            out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
            free(items);
        } else {
            out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        }
        free(path); free(visited);
    }

    graph_adj_free(a);
    return out;
}
