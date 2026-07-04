/* katzcentrality.c - KatzCentrality[g, alpha]: the Katz centrality of each
 * vertex with attenuation factor alpha (and base weight beta = 1). A vertex is
 * central if it is pointed to by central vertices, with the influence of a
 * length-k walk discounted by alpha^k:
 *
 *     x_i = alpha * sum_{j -> i} x_j + 1,   i.e.  (I - alpha A^T) x = 1.
 *
 * We assemble that linear system and solve it exactly through LinearSolve, so a
 * rational alpha yields an exact rational centrality vector. A^T uses in-edges
 * (predecessors), so for a directed graph a vertex's score depends on who points
 * at it; for an undirected graph in- and out-neighbourhoods coincide. alpha must
 * be a number (integer, rational, or real); the solve fails -- leaving the call
 * unevaluated -- when I - alpha A^T is singular (alpha at a reciprocal
 * eigenvalue). O(V^3) from the solve.
 *
 * Memory (SPEC section 4): returns a fresh vector; frees res. NULL on a
 * non-graph, a missing/non-numeric alpha, or a singular system.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include "eval.h"
#include <stdlib.h>
#include <string.h>

/* True if e is a plain numeric literal usable as alpha. */
static int is_number(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT)
        return 1;
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        strcmp(e->data.function.head->data.symbol, "Rational") == 0)
        return 1;
    return 0;
}

Expr* builtin_katz_centrality(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* alpha = res->data.function.args[1];
    if (!is_number(alpha)) return NULL;
    GraphAdj* a = graph_build_adj(res->data.function.args[0]);
    if (!a) return NULL;
    int n = a->n;
    if (n == 0) { graph_adj_free(a); return expr_new_function(expr_new_symbol(SYM_List), NULL, 0); }

    /* Predecessor adjacency: inadj[i][j] = 1 iff j -> i. */
    char* inadj = calloc((size_t)n * (size_t)n, 1);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < a->indeg[i]; j++) inadj[(size_t)i * n + a->in[i][j]] = 1;

    Expr** rows = calloc((size_t)n, sizeof(Expr*));
    Expr** bvec = calloc((size_t)n, sizeof(Expr*));
    for (int i = 0; i < n; i++) {
        Expr** row = calloc((size_t)n, sizeof(Expr*));
        for (int j = 0; j < n; j++) {
            if (i == j) {
                row[j] = expr_new_integer(1);                 /* diagonal (no self-loops) */
            } else if (inadj[(size_t)i * n + j]) {
                Expr* t[2] = { expr_new_integer(-1), expr_copy((Expr*)alpha) };
                row[j] = expr_new_function(expr_new_symbol(SYM_Times), t, 2);  /* -alpha */
            } else {
                row[j] = expr_new_integer(0);
            }
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), row, (size_t)n);
        free(row);
        bvec[i] = expr_new_integer(1);                        /* beta = 1 */
    }
    free(inadj);
    graph_adj_free(a);

    Expr* matrix = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    Expr* vec    = expr_new_function(expr_new_symbol(SYM_List), bvec, (size_t)n);
    free(rows); free(bvec);
    Expr* args[2] = { matrix, vec };
    Expr* solve = expr_new_function(expr_new_symbol("LinearSolve"), args, 2);
    return evaluate(solve);
}
