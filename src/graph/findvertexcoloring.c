/* findvertexcoloring.c - FindVertexColoring[g]: a proper vertex colouring of g
 * that uses as few colours as possible (the chromatic number), returned as a list
 * of colour indices 1..chi, one per vertex in vertex order.
 *
 * Tries k = 1, 2, ... colours and, for each, backtracks to colour the vertices so
 * that adjacent ones differ; the first feasible k is the chromatic number and its
 * colouring is returned. A symmetry cut (a vertex may open at most one brand-new
 * colour) keeps the search small. Edge direction is ignored. Adjacent vertices
 * always receive different colours; the number of distinct colours equals
 * ChromaticNumber[g].
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

typedef struct { const char* adj; int n; int k; int* color; } ColEnv;

static int color_from(ColEnv* e, int v, int maxUsed) {
    if (v == e->n) return 1;
    int limit = (maxUsed + 1 < e->k) ? maxUsed + 1 : e->k;
    for (int c = 1; c <= limit; c++) {
        int ok = 1;
        for (int u = 0; u < e->n; u++)
            if (e->adj[(size_t)v * e->n + u] && e->color[u] == c) { ok = 0; break; }
        if (!ok) continue;
        e->color[v] = c;
        if (color_from(e, v + 1, (c > maxUsed) ? c : maxUsed)) return 1;
        e->color[v] = 0;
    }
    return 0;
}

Expr* builtin_find_vertex_coloring(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    if (n == 0) return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);

    char* adj = calloc((size_t)n * (size_t)n, 1);
    for (size_t e = 0; e < edges->data.function.arg_count; e++) {
        const Expr* ed = edges->data.function.args[e];
        int ia = graph_vertex_index(verts, ed->data.function.args[0]);
        int ib = graph_vertex_index(verts, ed->data.function.args[1]);
        if (ia < 0 || ib < 0 || ia == ib) continue;
        adj[(size_t)ia * n + ib] = adj[(size_t)ib * n + ia] = 1;
    }

    int* color = calloc((size_t)n, sizeof(int));
    ColEnv e = { adj, n, 0, color };
    for (int k = 1; k <= n; k++) {
        e.k = k;
        for (int i = 0; i < n; i++) color[i] = 0;
        if (color_from(&e, 0, 0)) break;             /* minimal feasible k */
    }
    free(adj);

    Expr** items = calloc((size_t)n, sizeof(Expr*));
    for (int i = 0; i < n; i++) items[i] = expr_new_integer(color[i]);
    free(color);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
    free(items);
    return out;
}
