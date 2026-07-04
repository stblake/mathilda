/* findvertexcover.c - FindVertexCover[g]: a minimum vertex cover (a smallest set
 * of vertices touching every edge) as a vertex list, or {} when g has no edges.
 *
 * By the Gallai identity a minimum vertex cover is the complement of a maximum
 * independent set, so we find a maximum independent set exactly the same way as
 * FindIndependentVertexSet -- a maximum clique of the complement graph, via
 * branch-and-bound expansion over the non-adjacency -- and return the vertices
 * *outside* it. Every edge has at least one endpoint outside any independent
 * set, so the complement is a cover, and it is smallest because the independent
 * set is largest. Edge direction is ignored; deterministic. Exponential worst
 * case, tiny on real graphs.
 *
 * Unlike FindClique / FindIndependentVertexSet (which return a list containing
 * one set), FindVertexCover returns the cover as a flat vertex list, matching
 * Wolfram.
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

typedef struct {
    const char* nadj;   /* complement adjacency: nadj[i*n+j] = i,j non-adjacent */
    int n;
    int* best;
    int  bestlen;
} IndepEnv;

static void expand(IndepEnv* e, int* R, int rlen, int* P, int plen) {
    if (rlen > e->bestlen) { memcpy(e->best, R, (size_t)rlen * sizeof(int)); e->bestlen = rlen; }
    for (int i = 0; i < plen; i++) {
        if (rlen + (plen - i) <= e->bestlen) break;
        int v = P[i];
        R[rlen] = v;
        int* newP = (plen - i > 1) ? malloc((size_t)(plen - i - 1) * sizeof(int)) : NULL;
        int np = 0;
        for (int j = i + 1; j < plen; j++)
            if (e->nadj[(size_t)v * e->n + P[j]]) newP[np++] = P[j];
        expand(e, R, rlen + 1, newP, np);
        free(newP);
    }
}

Expr* builtin_find_vertex_cover(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    char* nadj = (n > 0) ? malloc((size_t)n * (size_t)n) : NULL;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            nadj[(size_t)i * n + j] = (i != j);
    for (size_t k = 0; k < ne; k++) {
        const Expr* ed = edges->data.function.args[k];
        int ia = graph_vertex_index(verts, ed->data.function.args[0]);
        int ib = graph_vertex_index(verts, ed->data.function.args[1]);
        if (ia < 0 || ib < 0 || ia == ib) continue;
        nadj[(size_t)ia * n + ib] = nadj[(size_t)ib * n + ia] = 0;
    }

    IndepEnv e = { nadj, n, (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL, 0 };
    if (n > 0) {
        int* R = malloc((size_t)n * sizeof(int));
        int* P = malloc((size_t)n * sizeof(int));
        for (int i = 0; i < n; i++) P[i] = i;
        expand(&e, R, 0, P, n);
        free(R); free(P);
    }
    free(nadj);

    /* cover = vertices not in the maximum independent set. */
    char* inset = (n > 0) ? calloc((size_t)n, 1) : NULL;
    for (int i = 0; i < e.bestlen; i++) inset[e.best[i]] = 1;
    int cover = n - e.bestlen;
    Expr** vs = (cover > 0) ? calloc((size_t)cover, sizeof(Expr*)) : NULL;
    int idx = 0;
    for (int v = 0; v < n; v++) if (!inset[v]) vs[idx++] = expr_copy(verts->data.function.args[v]);
    free(inset); free(e.best);

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), vs, (size_t)cover);
    free(vs);
    return out;
}
