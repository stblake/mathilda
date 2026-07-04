/* hamiltonianq.c - HamiltonianGraphQ[g]: True iff g has a Hamiltonian cycle (a
 * closed walk visiting every vertex exactly once). The predicate companion to
 * FindHamiltonianCycle, and the Hamiltonian counterpart of EulerianGraphQ.
 *
 * Depth-first backtracking over GraphAdj.out[] with visited-set pruning, rooted
 * at vertex 0 (WLOG -- a Hamiltonian cycle passes through every vertex). At full
 * depth the walk must close back to the start. Cheap necessary-condition prunes
 * short-circuit the trivially impossible: fewer than 3 vertices, or any vertex
 * missing an out- or in-neighbour. out[] encodes edge direction, so directed and
 * undirected graphs share one code path. O(V!) worst case, tiny in practice.
 *
 * Memory (SPEC section 4): returns True/False; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

static int ham_exists(const GraphAdj* a, int* path, char* visited, int depth, int start) {
    int last = path[depth - 1];
    if (depth == a->n) {
        for (int j = 0; j < a->outdeg[last]; j++)
            if (a->out[last][j] == start) return 1;
        return 0;
    }
    for (int j = 0; j < a->outdeg[last]; j++) {
        int u = a->out[last][j];
        if (!visited[u]) {
            visited[u] = 1;
            path[depth] = u;
            if (ham_exists(a, path, visited, depth + 1, start)) return 1;
            visited[u] = 0;
        }
    }
    return 0;
}

Expr* builtin_hamiltonian_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GraphAdj* a = graph_build_adj(res->data.function.args[0]);
    if (!a) return NULL;
    int n = a->n;

    int impossible = (n < 3);
    for (int i = 0; i < n && !impossible; i++)
        if (a->outdeg[i] == 0 || a->indeg[i] == 0) impossible = 1;

    int yes = 0;
    if (!impossible) {
        int*  path    = malloc((size_t)n * sizeof(int));
        char* visited = calloc((size_t)n, 1);
        path[0] = 0; visited[0] = 1;
        yes = ham_exists(a, path, visited, 1, 0);
        free(path); free(visited);
    }
    graph_adj_free(a);
    return expr_new_symbol(yes ? SYM_True : SYM_False);
}
