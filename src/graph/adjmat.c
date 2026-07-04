/* adjmat.c - AdjacencyMatrix[g]: dense 0/1 adjacency matrix.
 *
 * Returns an n x n dense List-of-Lists (n = |V|, in canonical vertex order),
 * consumable directly by Det, Tr, and Eigenvalues with no linalg changes.
 * A DirectedEdge[a,b] sets M[a][b] = 1; an UndirectedEdge sets both M[a][b] and
 * M[b][a], so undirected graphs yield a symmetric matrix. Entries are always
 * 0/1 since parallel edges are forbidden.
 *
 * Future hook: a WeightedAdjacencyMatrix would fill entries with edge weights
 * instead of 1 (Locked Decision 2); not implemented in the MVP.
 *
 * Memory (SPEC section 4): returns a freshly-allocated matrix; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_adjacency_matrix(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    size_t n = verts->data.function.arg_count;

    int* grid = (n > 0) ? calloc(n * n, sizeof(int)) : NULL;
    if (n > 0 && !grid) return NULL;

    for (size_t k = 0; k < edges->data.function.arg_count; k++) {
        const Expr* e = edges->data.function.args[k];
        const char* kind = graph_edge_kind(e);
        int ia = graph_vertex_index(verts, e->data.function.args[0]);
        int ib = graph_vertex_index(verts, e->data.function.args[1]);
        grid[(size_t)ia * n + (size_t)ib] = 1;
        if (kind == SYM_UndirectedEdge) grid[(size_t)ib * n + (size_t)ia] = 1;
    }

    Expr** rows = (n > 0) ? calloc(n, sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < n; i++) {
        Expr** row = calloc(n, sizeof(Expr*));
        for (size_t j = 0; j < n; j++)
            row[j] = expr_new_integer(grid[i * n + j]);
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), row, n);
        free(row);
    }
    Expr* mat = expr_new_function(expr_new_symbol(SYM_List), rows, n);
    free(rows);
    free(grid);
    return mat;
}
