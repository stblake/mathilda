/* findindependent.c - FindIndependentVertexSet[g]: a largest independent vertex
 * set (a set of pairwise non-adjacent vertices) as a list containing one set --
 * Wolfram's shape, e.g. {{2, 3, 4, 5}} for a star's leaves -- or {} when g has
 * no vertices.
 *
 * A maximum independent set of g is a maximum clique of the complement of g, so
 * we run the same branch-and-bound clique expansion (as in findclique.c) over
 * the complemented adjacency: candidates are the vertices *non-adjacent* to the
 * current set. Grow a set R over candidates P, keep the best, and prune the
 * branch when |R| + |remaining| <= |best|. Edge direction is ignored;
 * deterministic (candidates in vertex order, first maximum kept). Exponential
 * worst case, tiny on real graphs.
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

Expr* builtin_find_independent_vertex_set(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    /* Complement adjacency: start all-connected (off diagonal), clear real edges. */
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

    Expr* out;
    if (e.bestlen > 0) {
        char* inset = calloc((size_t)n, 1);
        for (int i = 0; i < e.bestlen; i++) inset[e.best[i]] = 1;
        Expr** vs = calloc((size_t)e.bestlen, sizeof(Expr*));
        int idx = 0;
        for (int v = 0; v < n; v++) if (inset[v]) vs[idx++] = expr_copy(verts->data.function.args[v]);
        free(inset);
        Expr* set = expr_new_function(expr_new_symbol(SYM_List), vs, (size_t)e.bestlen);
        free(vs);
        Expr* one[1] = { set };
        out = expr_new_function(expr_new_symbol(SYM_List), one, 1);
    } else {
        out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }
    free(e.best);
    return out;
}
