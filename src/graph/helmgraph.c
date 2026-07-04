/* helmgraph.c - HelmGraph[n]: the helm graph, a wheel (a hub joined to every
 * vertex of an n-cycle rim) with one pendant vertex attached to each rim vertex.
 *
 * Vertices: hub = 1, rim = 2..n+1 (an n-cycle), pendants = n+2..2n+1 (pendant r
 * attached to rim vertex r). Edges: the rim cycle, the n hub-rim spokes, and the
 * n pendants -- 2n+1 vertices and 3n edges. The hub has degree n, each rim vertex
 * degree 4 (two rim + spoke + pendant), each pendant degree 1. O(n).
 *
 * n >= 3 must be an integer. Memory (SPEC section 4): returns a freshly-built
 * Graph; frees res. NULL otherwise.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_helm_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* ne = res->data.function.args[0];
    if (ne->type != EXPR_INTEGER) return NULL;
    long n = (long)ne->data.integer;
    if (n < 3) return NULL;

    long V = 2 * n + 1;
    Expr** vc = malloc((size_t)V * sizeof(Expr*));
    for (long i = 0; i < V; i++) vc[i] = expr_new_integer(i + 1);

    Expr** ec = malloc((size_t)(3 * n) * sizeof(Expr*));
    size_t me = 0;
    #define ADD(a,b) do { Expr* _e[2] = { expr_new_integer(a), expr_new_integer(b) }; \
        ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), _e, 2); } while (0)
    for (long r = 0; r < n; r++) {
        long rim = 2 + r, rim1 = 2 + (r + 1) % n, pend = n + 2 + r;
        ADD(rim, rim1);      /* rim cycle */
        ADD(1, rim);         /* hub spoke */
        ADD(rim, pend);      /* pendant */
    }
    #undef ADD

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)V);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
