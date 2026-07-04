/* finddominatingset.c - FindDominatingSet[g]: a minimum dominating set of g -- a
 * smallest set S of vertices such that every vertex is in S or adjacent to a
 * vertex of S -- returned as a vertex list.
 *
 * Each vertex gets a closed-neighbourhood bitmask (itself plus its neighbours);
 * a set dominates iff the OR of its masks is all-ones. We search by increasing
 * size k = 1, 2, ... over the C(n,k) subsets and return the first dominating one,
 * so the result is a minimum dominating set. Edge direction is ignored. The
 * exponential subset search is bounded to modest n (DOM_MAX_N below, which also
 * keeps the vertex bitmasks within one 64-bit word); the call is left
 * unevaluated for larger graphs.
 *
 * A star is dominated by its centre; K_n by any single vertex; an edgeless graph
 * needs all its vertices. Memory (SPEC section 4): returns a fresh List; frees
 * res. NULL on a non-graph or an oversized graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

#define DOM_MAX_N 26

Expr* builtin_find_dominating_set(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    if (n == 0) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    if (n > DOM_MAX_N) return NULL;                       /* too large to search */

    /* Closed-neighbourhood masks (vertex i plus its undirected neighbours). */
    unsigned long long* mask = calloc((size_t)n, sizeof(unsigned long long));
    for (int i = 0; i < n; i++) mask[i] = 1ULL << i;
    for (size_t k = 0; k < edges->data.function.arg_count; k++) {
        const Expr* e = edges->data.function.args[k];
        int a = graph_vertex_index(verts, e->data.function.args[0]);
        int b = graph_vertex_index(verts, e->data.function.args[1]);
        if (a < 0 || b < 0 || a == b) continue;
        mask[a] |= 1ULL << b; mask[b] |= 1ULL << a;
    }
    unsigned long long full = (1ULL << n) - 1;   /* n <= DOM_MAX_N, so no overflow */

    /* Search subsets by increasing size; first dominating set is minimum. */
    int* idx = malloc((size_t)n * sizeof(int));
    int* best = NULL; int bestlen = -1;
    for (int k = 1; k <= n && bestlen < 0; k++) {
        for (int i = 0; i < k; i++) idx[i] = i;
        for (;;) {
            unsigned long long cover = 0;
            for (int i = 0; i < k; i++) cover |= mask[idx[i]];
            if (cover == full) { best = malloc((size_t)k * sizeof(int)); for (int i = 0; i < k; i++) best[i] = idx[i]; bestlen = k; break; }
            int p = k - 1;
            while (p >= 0 && idx[p] == n - k + p) p--;
            if (p < 0) break;
            idx[p]++;
            for (int i = p + 1; i < k; i++) idx[i] = idx[i - 1] + 1;
        }
    }
    free(idx); free(mask);

    Expr** items = (bestlen > 0) ? calloc((size_t)bestlen, sizeof(Expr*)) : NULL;
    for (int i = 0; i < bestlen; i++) items[i] = expr_copy(verts->data.function.args[best[i]]);
    free(best);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)(bestlen > 0 ? bestlen : 0));
    free(items);
    return out;
}
