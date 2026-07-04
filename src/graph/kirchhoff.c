/* kirchhoff.c - KirchhoffMatrix[g]: the graph Laplacian L = D - A.
 *
 * D is the diagonal matrix of vertex degrees (row sums of the adjacency matrix:
 * the degree for an undirected graph, the out-degree for a directed one) and A
 * is the 0/1 adjacency matrix. So L[i][i] = deg(i) and L[i][j] = -1 when there
 * is an edge i->j, else 0. Every row sums to 0.
 *
 * The Laplacian is the bridge to spectral graph theory: it is symmetric for an
 * undirected graph, the multiplicity of its zero eigenvalue is the number of
 * connected components, and (Matrix-Tree theorem) any cofactor equals the
 * number of spanning trees — all reachable via the existing dense linear
 * algebra (Eigenvalues, Det, ...) since the result is an ordinary matrix.
 *
 * O(V^2) to materialize the dense matrix. Memory (SPEC section 4): returns a
 * fresh List-of-Lists; frees res. NULL (unevaluated) on a non-graph argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_kirchhoff_matrix(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GraphAdj* a = graph_build_adj(res->data.function.args[0]);
    if (!a) return NULL;
    int n = a->n;

    /* Adjacency membership: adj[i*n+j] = 1 iff edge i->j. */
    char* adj = (n > 0) ? calloc((size_t)n * (size_t)n, 1) : NULL;
    if (n > 0 && !adj) { graph_adj_free(a); return NULL; }
    for (int i = 0; i < n; i++)
        for (int e = 0; e < a->outdeg[i]; e++)
            adj[(size_t)i * n + a->out[i][e]] = 1;

    Expr** rows = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    for (int i = 0; i < n; i++) {
        Expr** cells = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
        for (int j = 0; j < n; j++) {
            int v = (i == j) ? a->outdeg[i]
                             : (adj[(size_t)i * n + j] ? -1 : 0);
            cells[j] = expr_new_integer(v);
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)n);
        free(cells);
    }
    free(adj);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    graph_adj_free(a);
    return out;
}
