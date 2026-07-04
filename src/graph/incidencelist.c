/* incidencelist.c - IncidenceList[g, v]: the list of edges of g incident to
 * vertex v (v is one of the endpoints), in the graph's edge order. For a directed
 * graph this includes both the out-edges and in-edges at v. Complements
 * AdjacencyList, which returns the neighbouring vertices rather than the edges.
 *
 * O(E) with a linear endpoint comparison; edge kinds are preserved. An edge
 * appears once even if both endpoints equal v (impossible here -- no self-loops).
 * A vertex not in g yields {}.
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_incidence_list(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    const Expr* v = res->data.function.args[1];
    if (!graph_is_valid(g)) return NULL;

    const Expr* edges = g->data.function.args[1];
    size_t m = edges->data.function.arg_count;

    Expr** items = (m > 0) ? malloc(m * sizeof(Expr*)) : NULL;
    size_t cnt = 0;
    for (size_t k = 0; k < m; k++) {
        const Expr* e = edges->data.function.args[k];
        if (expr_eq(e->data.function.args[0], v) || expr_eq(e->data.function.args[1], v))
            items[cnt++] = expr_copy((Expr*)e);
    }

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, cnt);
    free(items);
    return out;
}
