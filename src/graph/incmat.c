/* incmat.c - IncidenceMatrix[g]: |V| x |E| incidence matrix.
 *
 * Column j corresponds to edge j (canonical order), row i to vertex i.
 *   - UndirectedEdge{a,b}: entries (a,j) and (b,j) are 1.
 *   - DirectedEdge[a,b]:   (a,j) = -1 (tail), (b,j) = 1 (head)  [oriented].
 *
 * Memory (SPEC section 4): returns a freshly-allocated matrix; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_incidence_matrix(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    size_t n = verts->data.function.arg_count;
    size_t m = edges->data.function.arg_count;

    int* grid = (n > 0 && m > 0) ? calloc(n * m, sizeof(int)) : NULL;
    if (n > 0 && m > 0 && !grid) return NULL;

    for (size_t j = 0; j < m; j++) {
        const Expr* e = edges->data.function.args[j];
        const char* kind = graph_edge_kind(e);
        int ia = graph_vertex_index(verts, e->data.function.args[0]);
        int ib = graph_vertex_index(verts, e->data.function.args[1]);
        if (kind == SYM_UndirectedEdge) {
            grid[(size_t)ia * m + j] = 1;
            grid[(size_t)ib * m + j] = 1;
        } else {
            grid[(size_t)ia * m + j] = -1;   /* tail */
            grid[(size_t)ib * m + j] = 1;    /* head */
        }
    }

    Expr** rows = (n > 0) ? calloc(n, sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < n; i++) {
        Expr** row = (m > 0) ? calloc(m, sizeof(Expr*)) : NULL;
        for (size_t j = 0; j < m; j++)
            row[j] = expr_new_integer(grid[i * m + j]);
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), row, m);
        free(row);
    }
    Expr* mat = expr_new_function(expr_new_symbol(SYM_List), rows, n);
    free(rows);
    free(grid);
    return mat;
}
