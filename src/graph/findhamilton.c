/* findhamilton.c - FindHamiltonianCycle[g]: a Hamiltonian cycle (a closed walk
 * visiting every vertex exactly once) as a vertex list {v0, v1, ..., v0}, or {}
 * if none exists. The constructive companion to the Eulerian-cycle finder.
 *
 * Depth-first backtracking with visited-set pruning, O(V!) worst case but tiny
 * in practice thanks to cheap necessary-condition prunes: any vertex missing an
 * out- or in-neighbour cannot lie on a cycle, so such graphs short-circuit to
 * {}. The search fixes the start at vertex 0 -- a Hamiltonian cycle, if it
 * exists, passes through every vertex, so this loses no generality and makes the
 * result deterministic. Successors come from GraphAdj.out[], which already
 * encodes edge direction (an undirected edge appears in both endpoints' out[]),
 * so the same code handles directed and undirected graphs.
 *
 * A Hamiltonian cycle needs at least 3 distinct vertices; with the graph's
 * simple-graph invariant (no self-loops, no parallel edges) every such cycle
 * automatically uses distinct edges, so no per-edge bookkeeping is required.
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Try to extend the partial path path[0..depth-1] into a full Hamiltonian cycle
 * rooted at `start`. Returns 1 and leaves path[0..n-1] filled on success. */
static int ham_extend(const GraphAdj* a, int* path, char* visited,
                      int depth, int start) {
    int last = path[depth - 1];
    if (depth == a->n) {                     /* all vertices used: close it? */
        for (int j = 0; j < a->outdeg[last]; j++)
            if (a->out[last][j] == start) return 1;
        return 0;
    }
    for (int j = 0; j < a->outdeg[last]; j++) {
        int u = a->out[last][j];
        if (!visited[u]) {
            visited[u] = 1;
            path[depth] = u;
            if (ham_extend(a, path, visited, depth + 1, start)) return 1;
            visited[u] = 0;
        }
    }
    return 0;
}

Expr* builtin_find_hamiltonian_cycle(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    int n = a->n;
    Expr* empty = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    /* Fewer than 3 vertices, or any vertex lacking an out/in neighbour, rules
     * out a cycle immediately. */
    int impossible = (n < 3);
    for (int i = 0; i < n && !impossible; i++)
        if (a->outdeg[i] == 0 || a->indeg[i] == 0) impossible = 1;
    if (impossible) { graph_adj_free(a); return empty; }

    int*  path    = malloc((size_t)n * sizeof(int));
    char* visited = calloc((size_t)n, 1);
    path[0] = 0; visited[0] = 1;

    Expr* out;
    if (ham_extend(a, path, visited, 1, 0)) {
        Expr** items = calloc((size_t)(n + 1), sizeof(Expr*));
        for (int i = 0; i < n; i++)
            items[i] = expr_copy((Expr*)a->verts->data.function.args[path[i]]);
        items[n] = expr_copy((Expr*)a->verts->data.function.args[path[0]]); /* close */
        out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)(n + 1));
        free(items);
        expr_free(empty);
    } else {
        out = empty;
    }

    free(path); free(visited);
    graph_adj_free(a);
    return out;
}
