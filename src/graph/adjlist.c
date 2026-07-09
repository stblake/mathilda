/* adjlist.c - AdjacencyList[g] and AdjacencyList[g, v].
 *
 * Neighbors of v: successors for directed edges (v -> u yields u), and both
 * endpoints for undirected edges. This matches the row convention of
 * AdjacencyMatrix (a 1 in row v, column u means v is adjacent to u). Neighbors
 * are returned in first-appearance order, de-duplicated.
 *
 * AdjacencyList[g]     -> {neighbors(v1), neighbors(v2), ...} in vertex order.
 * AdjacencyList[g, v]  -> neighbors(v).
 *
 * Memory (SPEC section 4): returns freshly-allocated lists; evaluator frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Build List of neighbors of `v` (deduplicated, first-appearance order). */
static Expr* neighbors_of(const Expr* g, const Expr* v) {
    const Expr* edges = g->data.function.args[1];
    size_t ne = edges->data.function.arg_count;
    Expr** nbr = (ne > 0) ? calloc(ne * 2, sizeof(Expr*)) : NULL;
    size_t n = 0;

    for (size_t i = 0; i < ne; i++) {
        const Expr* e = edges->data.function.args[i];
        const char* kind = graph_edge_kind(e);
        const Expr* a = e->data.function.args[0];
        const Expr* b = e->data.function.args[1];
        const Expr* add = NULL;
        if (kind == SYM_UndirectedEdge) {
            if (expr_eq(a, v)) add = b;
            else if (expr_eq(b, v)) add = a;
        } else {
            if (expr_eq(a, v)) add = b;   /* successor only */
        }
        if (add) {
            int seen = 0;
            for (size_t j = 0; j < n; j++)
                if (expr_eq(nbr[j], add)) { seen = 1; break; }
            if (!seen) nbr[n++] = expr_copy((Expr*)add);
        }
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), nbr, n);
    free(nbr);
    return list;
}

Expr* builtin_adjacency_list(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    const Expr* verts = g->data.function.args[0];

    if (argc == 2) {
        const Expr* v = res->data.function.args[1];
        if (graph_vertex_index(verts, v) < 0) return NULL;
        return neighbors_of(g, v);
    }

    size_t nv = verts->data.function.arg_count;
    Expr** rows = (nv > 0) ? calloc(nv, sizeof(Expr*)) : NULL;
    if (nv > 0 && !rows) return NULL;
    for (size_t i = 0; i < nv; i++)
        rows[i] = neighbors_of(g, verts->data.function.args[i]);
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), rows, nv);
    free(rows);
    return list;
}
