/* laddergraph.c - LadderGraph[n]: the ladder graph L_n, two paths on n vertices
 * ("rails") joined by n "rungs" -- the Cartesian product P_n [] P_2.
 *
 * Vertices are 1..2n: the top rail 1..n and the bottom rail n+1..2n. Edges are
 * the two rail paths (i,i+1) and (n+i,n+i+1) for i = 1..n-1, plus the rungs
 * (i, n+i) for i = 1..n -- 3n-2 edges in all. O(n). L_1 is a single edge, L_2 is
 * C_4 (a square).
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL unless
 * n >= 1 is an integer.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_ladder_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* ne = res->data.function.args[0];
    if (ne->type != EXPR_INTEGER) return NULL;
    long n = (long)ne->data.integer;
    if (n < 1) return NULL;

    long V = 2 * n;
    Expr** vc = malloc((size_t)V * sizeof(Expr*));
    for (long i = 0; i < V; i++) vc[i] = expr_new_integer(i + 1);

    size_t total = (size_t)(3 * n - 2);
    Expr** ec = malloc(total * sizeof(Expr*));
    size_t me = 0;
    for (long i = 1; i < n; i++) {                 /* top rail */
        Expr* a[2] = { expr_new_integer(i), expr_new_integer(i + 1) };
        ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), a, 2);
    }
    for (long i = 1; i < n; i++) {                 /* bottom rail */
        Expr* a[2] = { expr_new_integer(n + i), expr_new_integer(n + i + 1) };
        ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), a, 2);
    }
    for (long i = 1; i <= n; i++) {                /* rungs */
        Expr* a[2] = { expr_new_integer(i), expr_new_integer(n + i) };
        ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), a, 2);
    }

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)V);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
