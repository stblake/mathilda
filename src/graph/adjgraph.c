/* adjgraph.c - AdjacencyGraph[m] / AdjacencyGraph[{v1, ..., vn}, m]: build a
 * graph from a 0/1 adjacency matrix.
 *
 * The inverse of AdjacencyMatrix. Vertices are the integers 1..n, or the given
 * vertex list (which must have length n). A symmetric matrix yields an
 * undirected graph (one UndirectedEdge per i<j with m[i][j]=1); an asymmetric
 * matrix yields a directed graph (a DirectedEdge for each off-diagonal
 * m[i][j]=1). Diagonal entries (self-loops) are ignored. An entry other than
 * the exact integers 0 and 1 -- a multiplicity k >= 2 would be k parallel
 * edges, which Graph does not represent -- leaves the call unevaluated; the
 * entry is compared as the full int64 (or BigInt) value, never truncated.
 * Trailing options (DirectedEdges, EdgeWeight, EdgeCapacity) go through
 * graph_generator_with_options. The result is returned as a Graph[...]
 * expression and canonicalized/validated by the evaluator (builtin_graph).
 *
 * Memory (SPEC section 4): returns a freshly-allocated Graph; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>
#include "checked_alloc.h"

/* Matrix entry as 0 or 1, or -1 for anything else (another integer -- however
 * large -- a BigInt, or a non-integer). */
static int entry01(const Expr* row, size_t j) {
    const Expr* e = row->data.function.args[j];
    if (e->type != EXPR_INTEGER) return -1;
    if (e->data.integer == 0) return 0;
    if (e->data.integer == 1) return 1;
    return -1;
}

static Expr* adjacency_graph_core(Expr* res, int directed) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    const Expr* vin = (argc == 2) ? res->data.function.args[0] : NULL;
    const Expr* m = res->data.function.args[argc - 1];
    if (!graph_is_list(m)) return NULL;

    size_t n = m->data.function.arg_count;
    if (vin && (!graph_is_list(vin) || vin->data.function.arg_count != n))
        return NULL;
    /* Validate: n x n, entries 0/1. */
    for (size_t i = 0; i < n; i++) {
        const Expr* row = m->data.function.args[i];
        if (!graph_is_list(row) || row->data.function.arg_count != n) return NULL;
        for (size_t j = 0; j < n; j++)
            if (entry01(row, j) < 0) return NULL;
    }

    /* Symmetry test. */
    int symmetric = 1;
    for (size_t i = 0; i < n && symmetric; i++)
        for (size_t j = 0; j < n; j++)
            if (entry01(m->data.function.args[i], j)
                != entry01(m->data.function.args[j], i)) { symmetric = 0; break; }
    /* DirectedEdges -> True: a symmetric matrix is a pair of arcs per entry
     * pair, not one undirected edge; -> False cannot undo an asymmetry. */
    if (directed == 1) symmetric = 0;
    if (directed == 0 && !symmetric) return NULL;

    /* Vertex i (0-based) as an expression. */
#define VERT(i) (vin ? expr_copy(vin->data.function.args[(i)]) \
                     : expr_new_integer((int64_t)(i) + 1))
    Expr** verts = (n > 0) ? calloc(n, sizeof(Expr*)) : NULL;
    Expr** edges = (n > 0) ? calloc(n * n, sizeof(Expr*)) : NULL;
    if (n > 0 && (!verts || !edges)) {
        free(verts); free(edges);
        return mth_nomem_abort();
    }
    for (size_t i = 0; i < n; i++) verts[i] = VERT(i);
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, n);
    free(verts);

    /* Edges. Upper-bound count then trim via arg_count. */
    size_t ne = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i == j) continue;                    /* ignore self-loops */
            if (entry01(m->data.function.args[i], j) != 1) continue;
            if (symmetric && j < i) continue;        /* one undirected edge per pair */
            Expr* ea[2] = { VERT(i), VERT(j) };
            edges[ne++] = expr_new_function(
                expr_new_symbol(symmetric ? SYM_UndirectedEdge : SYM_DirectedEdge),
                ea, 2);
        }
    }
#undef VERT
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), edges, ne);
    free(edges);

    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}

Expr* builtin_adjacency_graph(Expr* res) {
    return graph_generator_with_options(res, 2, adjacency_graph_core);
}
