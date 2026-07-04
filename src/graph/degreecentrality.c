/* degreecentrality.c - DegreeCentrality[g]: for each vertex, the number of
 * edges incident to it -- the simplest centrality measure, completing the trio
 * with ClosenessCentrality and BetweennessCentrality.
 *
 * For an undirected graph this is the ordinary degree; for a directed graph it
 * is in-degree + out-degree. Both fall out of one rule -- every edge increments
 * a counter at each of its two endpoints -- since an undirected edge touches
 * each endpoint once (degree) and a directed edge a->b contributes one to a's
 * out-count and one to b's in-count. So no direction test is needed. O(V+E).
 *
 * Results are exact integers in vertex (canonical) order.
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_degree_centrality(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    long* deg = (n > 0) ? calloc((size_t)n, sizeof(long)) : NULL;
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        int ia = graph_vertex_index(verts, e->data.function.args[0]);
        int ib = graph_vertex_index(verts, e->data.function.args[1]);
        if (ia >= 0) deg[ia]++;
        if (ib >= 0) deg[ib]++;
    }

    Expr** items = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    for (int i = 0; i < n; i++) items[i] = expr_new_integer(deg[i]);
    free(deg);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
    free(items);
    return out;
}
