/* geargraph.c - GearGraph[n]: the gear (cogwheel) graph, a hub joined to every
 * other vertex of a 2n-cycle rim -- equivalently, a wheel with an extra vertex
 * inserted between each pair of adjacent rim vertices.
 *
 * Vertices: hub = 1 and a 2n-cycle rim 2..2n+1. Edges: the rim cycle (2n edges)
 * plus n spokes joining the hub to the even-indexed rim vertices. So 2n+1
 * vertices and 3n edges. It is bipartite (the even rim cycle 2-colours, and the
 * hub attaches only to one colour class). O(n).
 *
 * n >= 3 must be an integer. Memory (SPEC section 4): returns a freshly-built
 * Graph; frees res. NULL otherwise.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_gear_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* ne = res->data.function.args[0];
    if (ne->type != EXPR_INTEGER) return NULL;
    long n = (long)ne->data.integer;
    if (n < 3) return NULL;

    long rim = 2 * n, V = rim + 1;
    Expr** vc = malloc((size_t)V * sizeof(Expr*));
    for (long i = 0; i < V; i++) vc[i] = expr_new_integer(i + 1);

    Expr** ec = malloc((size_t)(3 * n) * sizeof(Expr*));
    size_t me = 0;
    #define ADD(a,b) do { Expr* _e[2] = { expr_new_integer(a), expr_new_integer(b) }; \
        ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), _e, 2); } while (0)
    for (long r = 0; r < rim; r++) {
        ADD(2 + r, 2 + (r + 1) % rim);         /* rim cycle */
        if (r % 2 == 0) ADD(1, 2 + r);         /* spoke to alternate rim vertices */
    }
    #undef ADD

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)V);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
