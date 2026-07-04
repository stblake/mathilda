/* degreesequence.c - DegreeSequence[g]: the vertex degrees sorted in
 * non-increasing order. For an undirected graph this is the ordinary degree
 * sequence; for a directed graph it is the total degree (in + out) per vertex,
 * sorted descending.
 *
 * Every edge increments a counter at each of its two endpoints (the same rule
 * as DegreeCentrality: an undirected edge touches each endpoint once = degree,
 * a directed edge a->b adds to a's out-count and b's in-count = in+out), then
 * the counts are sorted high-to-low. O(V + E + V log V). Exact integers.
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

static int cmp_desc(const void* a, const void* b) {
    long x = *(const long*)a, y = *(const long*)b;
    return (x < y) - (x > y);   /* descending */
}

Expr* builtin_degree_sequence(Expr* res) {
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
    if (n > 0) qsort(deg, (size_t)n, sizeof(long), cmp_desc);

    Expr** items = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    for (int i = 0; i < n; i++) items[i] = expr_new_integer(deg[i]);
    free(deg);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
    free(items);
    return out;
}
