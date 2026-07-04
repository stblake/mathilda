/* pagerank.c - PageRankCentrality[g]: the PageRank of each vertex -- the
 * stationary distribution of a random surfer that follows out-edges with
 * probability d = 17/20 and teleports uniformly with probability 1 - d.
 *
 * Instead of iterating to a floating-point fixed point, we solve the defining
 * linear system exactly and return rationals. PageRank satisfies
 *
 *     pi = (1-d)/n * 1  +  d * M pi,     i.e.  (I - d M) pi = (1-d)/n * 1,
 *
 * where M is the column-stochastic transition matrix: M[v][u] = 1/outdeg(u) if
 * u->v, and 1/n for a dangling u (no out-edges, teleport everywhere). Because M
 * is column-stochastic and d < 1, (I - dM) is invertible and the solution sums
 * to 1. We assemble the rational matrix and right-hand side and hand them to the
 * exact LinearSolve, so the whole result is a reduced-rational vector. Edge
 * direction is followed (out[] is symmetric for undirected graphs). O(V^3) from
 * the solve.
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include "eval.h"
#include <stdlib.h>

/* Exact num/den as a reduced value (integer or Rational). den > 0 required. */
static Expr* ratio(long num, long den) {
    Expr* pa[2] = { expr_new_integer(den), expr_new_integer(-1) };
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
    Expr* ta[2] = { expr_new_integer(num), inv };
    return evaluate(expr_new_function(expr_new_symbol(SYM_Times), ta, 2));
}

Expr* builtin_pagerank_centrality(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GraphAdj* a = graph_build_adj(res->data.function.args[0]);
    if (!a) return NULL;
    int n = a->n;
    if (n == 0) { graph_adj_free(a); return expr_new_function(expr_new_symbol(SYM_List), NULL, 0); }

    /* Successor adjacency as a boolean matrix. */
    char* out = calloc((size_t)n * (size_t)n, 1);
    for (int u = 0; u < n; u++)
        for (int j = 0; j < a->outdeg[u]; j++) out[(size_t)u * n + a->out[u][j]] = 1;

    /* Damping d = 17/20; teleport (1-d) = 3/20. */
    Expr** rows = calloc((size_t)n, sizeof(Expr*));
    Expr** bvec = calloc((size_t)n, sizeof(Expr*));
    for (int v = 0; v < n; v++) {
        Expr** row = calloc((size_t)n, sizeof(Expr*));
        for (int u = 0; u < n; u++) {
            long mnum, mden;
            if (a->outdeg[u] == 0) { mnum = 1; mden = n; }         /* dangling */
            else { mnum = out[(size_t)u * n + v] ? 1 : 0; mden = a->outdeg[u]; }
            /* A[v][u] = [v==u] - (17/20)*(mnum/mden)
             *         = ([v==u]*20*mden - 17*mnum) / (20*mden) */
            long num = (long)(v == u) * 20 * mden - 17 * mnum;
            row[u] = ratio(num, 20 * mden);
        }
        rows[v] = expr_new_function(expr_new_symbol(SYM_List), row, (size_t)n);
        free(row);
        bvec[v] = ratio(3, 20L * n);                                /* (1-d)/n */
    }
    free(out);
    graph_adj_free(a);

    Expr* matrix = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    Expr* vec    = expr_new_function(expr_new_symbol(SYM_List), bvec, (size_t)n);
    free(rows); free(bvec);
    Expr* args[2] = { matrix, vec };
    Expr* solve = expr_new_function(expr_new_symbol("LinearSolve"), args, 2);
    return evaluate(solve);
}
