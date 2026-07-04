/* edgebetweenness.c - EdgeBetweennessCentrality[g]: for each edge, the number of
 * shortest paths that run along it, summed over all vertex pairs -- the edge
 * analogue of BetweennessCentrality (and the score Girvan-Newman peels).
 *
 * The number of shortest s->t paths using directed edge a->b is sigma_sa *
 * sigma_bt whenever d(s,a) + 1 + d(b,t) = d(s,t); dividing by sigma_st gives the
 * fraction of shortest s-t paths through it. Summed over ordered pairs:
 *     eb(a->b) = sum_{s,t : d(s,a)+1+d(b,t)=d(s,t)} sigma_sa * sigma_bt / sigma_st.
 * An undirected edge carries paths both ways, so both a->b and b->a are summed
 * and the ordered total is halved (matching the undirected-pair convention of
 * BetweennessCentrality); a directed edge takes only its own orientation.
 *
 * All-pairs distances/path-counts come from one BFS per source (Brandes-style
 * counting); each term is num/sigma_st and the per-edge sum is reduced by the
 * evaluator, so results are exact rationals. O(E * V^2).
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include "eval.h"
#include <stdlib.h>

static void bfs_paths(const GraphAdj* a, int s, int* d, long* sig, int* q) {
    int n = a->n;
    for (int i = 0; i < n; i++) { d[i] = -1; sig[i] = 0; }
    int head = 0, tail = 0;
    d[s] = 0; sig[s] = 1; q[tail++] = s;
    while (head < tail) {
        int u = q[head++];
        for (int e = 0; e < a->outdeg[u]; e++) {
            int w = a->out[u][e];
            if (d[w] < 0) { d[w] = d[u] + 1; q[tail++] = w; }
            if (d[w] == d[u] + 1) sig[w] += sig[u];
        }
    }
}

/* Append the terms sigma_sa*sigma_bt/sigma_st for edge direction a->b. */
static void add_dir_terms(int n, const int* D, const long* SIG, int ai, int bi,
                          Expr** terms, size_t* k) {
    for (int s = 0; s < n; s++) {
        int dsa = D[(size_t)s * n + ai];
        if (dsa < 0) continue;
        long ssa = SIG[(size_t)s * n + ai];
        for (int t = 0; t < n; t++) {
            int dbt = D[(size_t)bi * n + t], dst = D[(size_t)s * n + t];
            if (dbt < 0 || dst < 0 || dsa + 1 + dbt != dst) continue;
            long num = ssa * SIG[(size_t)bi * n + t];
            long den = SIG[(size_t)s * n + t];
            if (num == 0 || den == 0) continue;
            Expr* pa[2] = { expr_new_integer(den), expr_new_integer(-1) };
            Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
            Expr* ta[2] = { expr_new_integer(num), inv };
            terms[(*k)++] = expr_new_function(expr_new_symbol(SYM_Times), ta, 2);
        }
    }
}

Expr* builtin_edge_betweenness_centrality(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t m = edges->data.function.arg_count;

    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;

    int*  D   = (n > 0) ? malloc((size_t)n * n * sizeof(int))  : NULL;
    long* SIG = (n > 0) ? malloc((size_t)n * n * sizeof(long)) : NULL;
    int*  q   = (n > 0) ? malloc((size_t)n * sizeof(int))      : NULL;
    for (int s = 0; s < n; s++) bfs_paths(a, s, D + (size_t)s * n, SIG + (size_t)s * n, q);
    free(q);
    graph_adj_free(a);

    Expr** items = (m > 0) ? calloc(m, sizeof(Expr*)) : NULL;
    for (size_t e = 0; e < m; e++) {
        const Expr* ed = edges->data.function.args[e];
        int undirected = (graph_edge_kind(ed) == SYM_UndirectedEdge);
        int ai = graph_vertex_index(verts, ed->data.function.args[0]);
        int bi = graph_vertex_index(verts, ed->data.function.args[1]);
        Expr** terms = (n > 0) ? calloc((size_t)2 * n * n, sizeof(Expr*)) : NULL;
        size_t k = 0;
        add_dir_terms(n, D, SIG, ai, bi, terms, &k);
        if (undirected) add_dir_terms(n, D, SIG, bi, ai, terms, &k);
        Expr* sum = (k == 0) ? expr_new_integer(0)
                             : expr_new_function(expr_new_symbol(SYM_Plus), terms, k);
        free(terms);
        if (undirected && k > 0) {                 /* halve the ordered-pair sum */
            Expr* ha[2] = { expr_new_integer(2), expr_new_integer(-1) };
            Expr* half = expr_new_function(expr_new_symbol(SYM_Power), ha, 2);
            Expr* ma[2] = { sum, half };
            sum = expr_new_function(expr_new_symbol(SYM_Times), ma, 2);
        }
        items[e] = evaluate(sum);
    }
    free(D); free(SIG);

    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, m);
    free(items);
    return out;
}
