/* findclique.c - FindClique[g]: a largest clique (a set of pairwise-adjacent
 * vertices) as a list containing one clique -- Wolfram's shape, e.g.
 * {{1, 2, 3, 4}} for K4 -- or {} for a graph with no vertices.
 *
 * Recursive clique expansion (a max-clique specialisation of Bron-Kerbosch): grow
 * a clique R over a candidate set P of vertices adjacent to all of R, considering
 * only candidates after the last-added one so each clique is reached once. The
 * best clique seen is kept; the branch `|R| + |remaining candidates| <= |best|`
 * is pruned, which keeps the exponential worst case tiny on real graphs. Edge
 * direction is ignored (the clique is defined on the underlying undirected
 * graph). Deterministic: candidates in vertex order, the first maximum kept.
 *
 * Memory (SPEC section 4): returns a fresh List of one List; frees res. NULL on
 * a non-graph argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

typedef struct {
    const char* adj;
    int n;
    int* best;
    int  bestlen;
} CliqueEnv;

static void expand(CliqueEnv* e, int* R, int rlen, int* P, int plen) {
    if (rlen > e->bestlen) { memcpy(e->best, R, (size_t)rlen * sizeof(int)); e->bestlen = rlen; }
    for (int i = 0; i < plen; i++) {
        if (rlen + (plen - i) <= e->bestlen) break;      /* can't beat best */
        int v = P[i];
        R[rlen] = v;
        int* newP = (plen - i > 1) ? malloc((size_t)(plen - i - 1) * sizeof(int)) : NULL;
        int np = 0;
        for (int j = i + 1; j < plen; j++)
            if (e->adj[(size_t)v * e->n + P[j]]) newP[np++] = P[j];
        expand(e, R, rlen + 1, newP, np);
        free(newP);
    }
}

Expr* builtin_find_clique(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    char* adj = (n > 0) ? calloc((size_t)n * (size_t)n, 1) : NULL;
    for (size_t k = 0; k < ne; k++) {
        const Expr* ed = edges->data.function.args[k];
        int ia = graph_vertex_index(verts, ed->data.function.args[0]);
        int ib = graph_vertex_index(verts, ed->data.function.args[1]);
        if (ia < 0 || ib < 0 || ia == ib) continue;
        adj[(size_t)ia * n + ib] = adj[(size_t)ib * n + ia] = 1;
    }

    CliqueEnv e = { adj, n, (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL, 0 };
    if (n > 0) {
        int* R = malloc((size_t)n * sizeof(int));
        int* P = malloc((size_t)n * sizeof(int));
        for (int i = 0; i < n; i++) P[i] = i;
        expand(&e, R, 0, P, n);
        free(R); free(P);
    }
    free(adj);

    Expr* out;
    if (e.bestlen > 0) {
        /* best[] holds indices in discovery order; emit in canonical order. */
        char* inset = calloc((size_t)n, 1);
        for (int i = 0; i < e.bestlen; i++) inset[e.best[i]] = 1;
        Expr** vs = calloc((size_t)e.bestlen, sizeof(Expr*));
        int idx = 0;
        for (int v = 0; v < n; v++) if (inset[v]) vs[idx++] = expr_copy(verts->data.function.args[v]);
        free(inset);
        Expr* clique = expr_new_function(expr_new_symbol(SYM_List), vs, (size_t)e.bestlen);
        free(vs);
        Expr* one[1] = { clique };
        out = expr_new_function(expr_new_symbol(SYM_List), one, 1);
    } else {
        out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }
    free(e.best);
    return out;
}
