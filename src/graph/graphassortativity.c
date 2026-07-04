/* graphassortativity.c - GraphAssortativity[g]: the degree assortativity
 * coefficient of g -- the Pearson correlation between the degrees of the two
 * endpoints of an edge, in [-1, 1]. It is +1 when high-degree vertices attach to
 * high-degree ones and -1 when they attach to low-degree ones.
 *
 * Newman's edge-list form reduces to a ratio of integers, so the result is an
 * exact rational. Over the M undirected edges, with j,k the endpoint degrees, let
 *   A   = sum j*k,   S1 = sum (j+k),   Sq = sum (j^2 + k^2).
 * Then r = (4*M*A - S1^2) / (2*M*Sq - S1^2). The denominator is the degree
 * variance times a positive factor; it is 0 exactly when g is regular (all
 * endpoint degrees equal), in which case the coefficient is Indeterminate.
 *
 * Edge direction is ignored (computed on the underlying undirected graph). A star
 * is perfectly disassortative (-1). Memory (SPEC section 4): returns a fresh
 * number; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include "eval.h"
#include <stdlib.h>

static Expr* ratio(long num, long den) {
    Expr* pa[2] = { expr_new_integer(den), expr_new_integer(-1) };
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
    Expr* ta[2] = { expr_new_integer(num), inv };
    return evaluate(expr_new_function(expr_new_symbol(SYM_Times), ta, 2));
}

Expr* builtin_graph_assortativity(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    /* Undirected adjacency (dedup) + degrees. */
    char* adj = (n > 0) ? calloc((size_t)n * (size_t)n, 1) : NULL;
    int*  deg = (n > 0) ? calloc((size_t)n, sizeof(int)) : NULL;
    for (size_t e = 0; e < ne; e++) {
        const Expr* ed = edges->data.function.args[e];
        int a = graph_vertex_index(verts, ed->data.function.args[0]);
        int b = graph_vertex_index(verts, ed->data.function.args[1]);
        if (a < 0 || b < 0 || a == b) continue;
        if (!adj[(size_t)a * n + b]) {
            adj[(size_t)a * n + b] = adj[(size_t)b * n + a] = 1;
            deg[a]++; deg[b]++;
        }
    }

    long M = 0, A = 0, S1 = 0, Sq = 0;
    for (int i = 0; i < n; i++)
        for (int k = i + 1; k < n; k++)
            if (adj[(size_t)i * n + k]) {
                long j = deg[i], l = deg[k];
                M++; A += j * l; S1 += j + l; Sq += j * j + l * l;
            }
    free(adj); free(deg);

    long den = 2 * M * Sq - S1 * S1;
    if (den == 0) return expr_new_symbol("Indeterminate");   /* regular / edgeless */
    return ratio(4 * M * A - S1 * S1, den);
}
