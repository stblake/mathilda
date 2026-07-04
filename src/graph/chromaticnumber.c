/* chromaticnumber.c - ChromaticNumber[g]: the least number of colours needed to
 * properly colour g (adjacent vertices differ), i.e. the smallest k for which a
 * proper k-colouring exists.
 *
 * Rather than reading it off the chromatic polynomial (exponential in edges),
 * this tries k = 1, 2, ... and tests k-colourability directly by backtracking:
 * colour vertices in order, giving each a colour in 1..k that clashes with no
 * already-coloured neighbour. Two symmetry cuts keep it fast: a vertex may open
 * at most one brand-new colour (colour <= maxUsed + 1), and the first k that
 * succeeds is the chromatic number, so the search stops immediately. Edge
 * direction is ignored (colouring is an undirected notion).
 *
 * Conventions: an empty graph (no vertices) has chromatic number 0; a graph with
 * vertices but no edges needs 1 colour; a bipartite graph 2; K_n needs n.
 *
 * Memory (SPEC section 4): returns a fresh integer; frees res. NULL on a
 * non-graph argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

typedef struct { const char* adj; int n; int k; int* color; } ColorEnv;

/* Try to colour vertices v..n-1; maxUsed = highest colour used in 0..v-1. */
static int color_from(ColorEnv* e, int v, int maxUsed) {
    if (v == e->n) return 1;
    int limit = (maxUsed + 1 < e->k) ? maxUsed + 1 : e->k;   /* symmetry cut */
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

Expr* builtin_chromatic_number(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;
    if (n == 0) return expr_new_integer(0);

    char* adj = calloc((size_t)n * (size_t)n, 1);
    int has_edge = 0;
    for (size_t k = 0; k < ne; k++) {
        const Expr* ed = edges->data.function.args[k];
        int ia = graph_vertex_index(verts, ed->data.function.args[0]);
        int ib = graph_vertex_index(verts, ed->data.function.args[1]);
        if (ia < 0 || ib < 0 || ia == ib) continue;
        adj[(size_t)ia * n + ib] = adj[(size_t)ib * n + ia] = 1;
        has_edge = 1;
    }
    if (!has_edge) { free(adj); return expr_new_integer(1); }

    int* color = calloc((size_t)n, sizeof(int));
    ColorEnv e = { adj, n, 0, color };
    int chi = n;                       /* K_n upper bound; loop finds the least */
    for (int k = 2; k <= n; k++) {
        e.k = k;
        for (int i = 0; i < n; i++) color[i] = 0;
        if (color_from(&e, 0, 0)) { chi = k; break; }
    }
    free(color); free(adj);
    return expr_new_integer(chi);
}
