/* adjgraph.c - AdjacencyGraph[m]: build a graph from a 0/1 adjacency matrix.
 *
 * The inverse of AdjacencyMatrix. Vertices are the integers 1..n. A symmetric
 * matrix yields an undirected graph (one UndirectedEdge per i<j with m[i][j]=1);
 * an asymmetric matrix yields a directed graph (a DirectedEdge for each off-
 * diagonal m[i][j]=1). Diagonal entries (self-loops) are ignored. The result is
 * returned as a Graph[...] expression and canonicalized/validated by the
 * evaluator (builtin_graph).
 *
 * Round-trips with AdjacencyMatrix when the source graph's vertices are 1..n.
 *
 * Memory (SPEC section 4): returns a freshly-allocated Graph; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Integer value of matrix entry, or -1 if not a plain integer. */
static int entry_int(const Expr* row, size_t j) {
    const Expr* e = row->data.function.args[j];
    if (e->type != EXPR_INTEGER) return -1;
    return (int)e->data.integer;
}

Expr* builtin_adjacency_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* m = res->data.function.args[0];
    if (!graph_is_list(m)) return NULL;

    size_t n = m->data.function.arg_count;
    /* Validate: n x n, entries 0/1. */
    for (size_t i = 0; i < n; i++) {
        const Expr* row = m->data.function.args[i];
        if (!graph_is_list(row) || row->data.function.arg_count != n) return NULL;
        for (size_t j = 0; j < n; j++) {
            int v = entry_int(row, j);
            if (v != 0 && v != 1) return NULL;
        }
    }

    /* Symmetry test. */
    int symmetric = 1;
    for (size_t i = 0; i < n && symmetric; i++)
        for (size_t j = 0; j < n; j++)
            if (entry_int(m->data.function.args[i], j)
                != entry_int(m->data.function.args[j], i)) { symmetric = 0; break; }

    /* Vertices 1..n. */
    Expr** verts = (n > 0) ? calloc(n, sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < n; i++) verts[i] = expr_new_integer((int64_t)i + 1);
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, n);
    free(verts);

    /* Edges. Upper-bound count then trim via arg_count. */
    Expr** edges = (n > 0) ? calloc(n * n, sizeof(Expr*)) : NULL;
    size_t ne = 0;
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i == j) continue;                    /* ignore self-loops */
            if (entry_int(m->data.function.args[i], j) != 1) continue;
            if (symmetric && j < i) continue;        /* one undirected edge per pair */
            Expr* ea[2] = { expr_new_integer((int64_t)i + 1),
                            expr_new_integer((int64_t)j + 1) };
            edges[ne++] = expr_new_function(
                expr_new_symbol(symmetric ? SYM_UndirectedEdge : SYM_DirectedEdge),
                ea, 2);
        }
    }
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), edges, ne);
    free(edges);

    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
