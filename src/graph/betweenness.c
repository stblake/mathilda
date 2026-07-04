/* betweenness.c - BetweennessCentrality[g].
 *
 * For a vertex v, the betweenness is the sum over ordered pairs (s,t), both != v,
 * of the fraction of shortest s-t paths that pass through v:
 *     c(v) = sum_{s != v != t} sigma_sv * sigma_vt / sigma_st ,   over pairs with
 *            d(s,v) + d(v,t) = d(s,t),
 * where sigma_ab counts shortest a->b paths. Undirected graphs count each
 * unordered pair once (the ordered sum is halved), matching the Wolfram Language;
 * directed graphs keep the ordered sum. Distances follow edge direction.
 *
 * All-pairs distances and path counts come from one BFS per source (Brandes-style
 * counting), then an O(V^3) accumulation. Values are exact: each term is
 * num/sigma_st and the per-vertex sum is reduced by the evaluator, so a graph
 * with tied shortest paths yields clean fractions (every C4 vertex is 1/2).
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include "eval.h"
#include <stdlib.h>

/* BFS from s over out[]: fills d[] (-1 unreachable) and sig[] (# shortest paths). */
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

Expr* builtin_betweenness_centrality(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    const Expr* edges = g->data.function.args[1];
    size_t ne = edges->data.function.arg_count;
    int directed = (ne > 0);
    for (size_t i = 0; i < ne; i++)
        if (graph_edge_kind(edges->data.function.args[i]) != SYM_DirectedEdge) { directed = 0; break; }

    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    int n = a->n;

    int*  D   = (n > 0) ? malloc((size_t)n * n * sizeof(int))  : NULL;
    long* SIG = (n > 0) ? malloc((size_t)n * n * sizeof(long)) : NULL;
    int*  q   = (n > 0) ? malloc((size_t)n * sizeof(int))      : NULL;
    if (n > 0 && (!D || !SIG || !q)) { free(D); free(SIG); free(q); graph_adj_free(a); return NULL; }
    for (int s = 0; s < n; s++) bfs_paths(a, s, D + (size_t)s * n, SIG + (size_t)s * n, q);
    free(q);

    Expr** items = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    for (int v = 0; v < n; v++) {
        /* Collect the nonzero terms sigma_sv * sigma_vt / sigma_st. */
        Expr** terms = (n > 1) ? calloc((size_t)n * n, sizeof(Expr*)) : NULL;
        size_t k = 0;
        for (int s = 0; s < n; s++) {
            if (s == v) continue;
            long dsv = D[(size_t)s * n + v], ssv = SIG[(size_t)s * n + v];
            if (dsv < 0) continue;
            for (int t = 0; t < n; t++) {
                if (t == v || t == s) continue;
                long dvt = D[(size_t)v * n + t], dst = D[(size_t)s * n + t];
                if (dvt < 0 || dst < 0 || dsv + dvt != dst) continue;
                long num = ssv * SIG[(size_t)v * n + t];   /* sigma_sv * sigma_vt */
                long den = SIG[(size_t)s * n + t];          /* sigma_st           */
                if (num == 0 || den == 0) continue;
                Expr* pa[2] = { expr_new_integer(den), expr_new_integer(-1) };
                Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
                Expr* ta[2] = { expr_new_integer(num), inv };
                terms[k++] = expr_new_function(expr_new_symbol(SYM_Times), ta, 2);
            }
        }
        Expr* sum = (k == 0) ? expr_new_integer(0)
                             : expr_new_function(expr_new_symbol(SYM_Plus), terms, k);
        free(terms);
        if (!directed && k > 0) {   /* halve the ordered-pair sum */
            Expr* ha[2] = { expr_new_integer(2), expr_new_integer(-1) };
            Expr* half = expr_new_function(expr_new_symbol(SYM_Power), ha, 2);
            Expr* ma[2] = { sum, half };
            sum = expr_new_function(expr_new_symbol(SYM_Times), ma, 2);
        }
        items[v] = evaluate(sum);
    }
    free(D); free(SIG);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)n);
    free(items);
    graph_adj_free(a);
    return out;
}
