/* wtadjmat.c - WeightedAdjacencyMatrix[g]: dense adjacency matrix filled with
 * per-edge weights instead of a literal 1.
 *
 * Same algorithm as AdjacencyMatrix (adjmat.c): a DirectedEdge[a,b] sets
 * M[a][b] = weight(a,b); an UndirectedEdge sets both M[a][b] and M[b][a]. Any
 * entry with no edge is 0. For a graph with no EdgeWeight, every weight
 * defaults to 1 (graph_resolve_edge_weights), so
 * WeightedAdjacencyMatrix[g] == AdjacencyMatrix[g] exactly for an unweighted g.
 *
 * Memory (SPEC section 4): returns a freshly-allocated matrix; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_weighted_adjacency_matrix(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    size_t n = verts->data.function.arg_count;

    Expr* weights = graph_resolve_edge_weights(g);
    if (!weights) return NULL;

    Expr** grid = (n > 0) ? calloc(n * n, sizeof(Expr*)) : NULL;
    if (n > 0 && !grid) { expr_free(weights); return NULL; }

    GraphVIdx* ix = graph_vidx_new(n);
    if (!ix) { free(grid); expr_free(weights); return NULL; }
    for (size_t i = 0; i < n; i++)
        graph_vidx_put(ix, verts->data.function.args[i], (int)i);

    for (size_t k = 0; k < edges->data.function.arg_count; k++) {
        const Expr* e = edges->data.function.args[k];
        const char* kind = graph_edge_kind(e);
        int ia = graph_vidx_get(ix, e->data.function.args[0]);
        int ib = graph_vidx_get(ix, e->data.function.args[1]);
        if (ia < 0 || ib < 0) continue;              /* validated: cannot happen */
        Expr* w = weights->data.function.args[k];
        grid[(size_t)ia * n + (size_t)ib] = expr_copy(w);
        if (kind == SYM_UndirectedEdge) grid[(size_t)ib * n + (size_t)ia] = expr_copy(w);
    }
    graph_vidx_free(ix);
    expr_free(weights);

    Expr** rows = (n > 0) ? calloc(n, sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < n; i++) {
        Expr** row = calloc(n, sizeof(Expr*));
        for (size_t j = 0; j < n; j++)
            row[j] = grid[i * n + j] ? grid[i * n + j] : expr_new_integer(0);
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), row, n);
        free(row);
    }
    Expr* mat = expr_new_function(expr_new_symbol(SYM_List), rows, n);
    free(rows);
    free(grid);
    return mat;
}
