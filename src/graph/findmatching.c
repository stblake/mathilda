/* findmatching.c - FindIndependentEdgeSet[g]: a maximum matching of g -- a
 * largest set of edges no two of which share a vertex -- returned as a list of
 * edges.
 *
 * Depth-first branch-and-bound over the edges: each edge is either taken (if both
 * endpoints are still free) or skipped, tracking the best matching found and
 * pruning when the current size plus the untried edges cannot beat it. Edge
 * direction is irrelevant to independence (two edges conflict iff they share an
 * endpoint). Exponential worst case, tiny on real graphs; deterministic (edges in
 * order, first maximum kept).
 *
 * K_{2k} and even cycles/paths yield perfect/near-perfect matchings; a star has a
 * maximum matching of size 1. Memory (SPEC section 4): returns a fresh List;
 * frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

typedef struct {
    const int* ea; const int* eb;    /* edge endpoint vertex indices */
    int m;
    char* used;                       /* vertex occupied by the current matching */
    int* cur; int curlen;             /* current matching (edge indices)          */
    int* best; int bestlen;
} MatchEnv;

static void expand(MatchEnv* e, int i, int count) {
    if (count > e->bestlen) { memcpy(e->best, e->cur, (size_t)count * sizeof(int)); e->bestlen = count; }
    if (i == e->m || count + (e->m - i) <= e->bestlen) return;   /* done / bound */
    int u = e->ea[i], v = e->eb[i];
    if (!e->used[u] && !e->used[v]) {                            /* take edge i */
        e->used[u] = e->used[v] = 1;
        e->cur[count] = i;
        expand(e, i + 1, count + 1);
        e->used[u] = e->used[v] = 0;
    }
    expand(e, i + 1, count);                                     /* skip edge i */
}

Expr* builtin_find_independent_edge_set(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    int m = (int)edges->data.function.arg_count;

    int* ea = (m > 0) ? malloc((size_t)m * sizeof(int)) : NULL;
    int* eb = (m > 0) ? malloc((size_t)m * sizeof(int)) : NULL;
    for (int k = 0; k < m; k++) {
        const Expr* e = edges->data.function.args[k];
        ea[k] = graph_vertex_index(verts, e->data.function.args[0]);
        eb[k] = graph_vertex_index(verts, e->data.function.args[1]);
    }

    MatchEnv env;
    env.ea = ea; env.eb = eb; env.m = m;
    env.used = (n > 0) ? calloc((size_t)n, 1) : NULL;
    env.cur  = (m > 0) ? malloc((size_t)m * sizeof(int)) : NULL;
    env.best = (m > 0) ? malloc((size_t)m * sizeof(int)) : NULL;
    env.curlen = env.bestlen = 0;
    if (m > 0) expand(&env, 0, 0);

    Expr** items = (env.bestlen > 0) ? calloc((size_t)env.bestlen, sizeof(Expr*)) : NULL;
    for (int i = 0; i < env.bestlen; i++)
        items[i] = expr_copy(edges->data.function.args[env.best[i]]);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)env.bestlen);

    free(ea); free(eb); free(env.used); free(env.cur); free(env.best); free(items);
    return out;
}
