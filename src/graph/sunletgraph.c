/* sunletgraph.c - SunletGraph[n]: the n-sunlet graph, a cycle C_n with one
 * pendant vertex attached to each cycle vertex ("rays of the sun").
 *
 * Vertices are the cycle 1..n and the pendants n+1..2n; edges are the cycle
 * i ~ i+1 (mod n) and the pendant edges i ~ n+i. So 2n vertices and 2n edges:
 * each cycle vertex has degree 3, each pendant degree 1. O(n). Also called the
 * sun/corona C_n o K_1.
 *
 * n >= 3 must be an integer. Memory (SPEC section 4): returns a freshly-built
 * Graph; frees res. NULL otherwise.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_sunlet_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* ne = res->data.function.args[0];
    if (ne->type != EXPR_INTEGER) return NULL;
    long n = (long)ne->data.integer;
    if (n < 3) return NULL;

    long V = 2 * n;
    Expr** vc = malloc((size_t)V * sizeof(Expr*));
    for (long i = 0; i < V; i++) vc[i] = expr_new_integer(i + 1);

    Expr** ec = malloc((size_t)(2 * n) * sizeof(Expr*));
    size_t me = 0;
    #define ADD(a,b) do { Expr* _e[2] = { expr_new_integer(a), expr_new_integer(b) }; \
        ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), _e, 2); } while (0)
    for (long i = 0; i < n; i++) {
        ADD(i + 1, (i + 1) % n + 1);   /* cycle edge */
        ADD(i + 1, n + i + 1);         /* pendant edge */
    }
    #undef ADD

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)V);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
