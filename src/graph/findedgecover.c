/* findedgecover.c - FindEdgeCover[g]: a minimum edge cover -- a smallest set of
 * edges such that every vertex is incident to at least one of them -- returned as
 * a list of edges.
 *
 * By Gallai's construction a minimum edge cover has size n - (maximum matching),
 * and is obtained by taking a maximum matching M and adding, for each vertex not
 * covered by M, one edge incident to it. (Unmatched vertices are pairwise
 * non-adjacent -- else M would not be maximum -- so each added edge covers
 * exactly one new vertex, giving the minimum.) A maximum matching is found by the
 * same branch-and-bound as FindIndependentEdgeSet.
 *
 * An edge cover exists iff g has no isolated vertex; when one does, {} is returned
 * (no cover). Edge direction is irrelevant. Memory (SPEC section 4): returns a
 * fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

typedef struct {
    const int* ea; const int* eb; int m;
    char* used; int* cur; int* best; int bestlen;
} MEnv;

static void expand(MEnv* e, int i, int count) {
    if (count > e->bestlen) { memcpy(e->best, e->cur, (size_t)count * sizeof(int)); e->bestlen = count; }
    if (i == e->m || count + (e->m - i) <= e->bestlen) return;
    int u = e->ea[i], v = e->eb[i];
    if (!e->used[u] && !e->used[v]) {
        e->used[u] = e->used[v] = 1; e->cur[count] = i;
        expand(e, i + 1, count + 1);
        e->used[u] = e->used[v] = 0;
    }
    expand(e, i + 1, count);
}

Expr* builtin_find_edge_cover(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    int m = (int)edges->data.function.arg_count;
    Expr* empty = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    if (n == 0) return empty;

    int* ea = malloc((size_t)(m > 0 ? m : 1) * sizeof(int));
    int* eb = malloc((size_t)(m > 0 ? m : 1) * sizeof(int));
    int* deg = calloc((size_t)n, sizeof(int));
    for (int k = 0; k < m; k++) {
        const Expr* e = edges->data.function.args[k];
        ea[k] = graph_vertex_index(verts, e->data.function.args[0]);
        eb[k] = graph_vertex_index(verts, e->data.function.args[1]);
        if (ea[k] >= 0) deg[ea[k]]++;
        if (eb[k] >= 0) deg[eb[k]]++;
    }
    for (int i = 0; i < n; i++) if (deg[i] == 0) {   /* isolated -> no cover */
        free(ea); free(eb); free(deg); return empty;
    }
    free(deg);

    /* Maximum matching. */
    MEnv env; env.ea = ea; env.eb = eb; env.m = m;
    env.used = calloc((size_t)n, 1);
    env.cur  = malloc((size_t)(m > 0 ? m : 1) * sizeof(int));
    env.best = malloc((size_t)(m > 0 ? m : 1) * sizeof(int));
    env.bestlen = 0;
    if (m > 0) expand(&env, 0, 0);

    /* Cover = matching edges, then one incident edge per uncovered vertex. */
    char* chosen  = calloc((size_t)(m > 0 ? m : 1), 1);
    char* covered = calloc((size_t)n, 1);
    for (int i = 0; i < env.bestlen; i++) {
        int k = env.best[i]; chosen[k] = 1;
        covered[ea[k]] = covered[eb[k]] = 1;
    }
    for (int v = 0; v < n; v++) {
        if (covered[v]) continue;
        for (int k = 0; k < m; k++)
            if (ea[k] == v || eb[k] == v) { chosen[k] = 1; covered[ea[k]] = covered[eb[k]] = 1; break; }
    }

    size_t cnt = 0;
    for (int k = 0; k < m; k++) cnt += chosen[k];
    Expr** items = (cnt > 0) ? calloc(cnt, sizeof(Expr*)) : NULL;
    size_t idx = 0;
    for (int k = 0; k < m; k++) if (chosen[k]) items[idx++] = expr_copy(edges->data.function.args[k]);
    free(ea); free(eb); free(env.used); free(env.cur); free(env.best); free(chosen); free(covered);
    expr_free(empty);

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, cnt);
    free(items);
    return out;
}
