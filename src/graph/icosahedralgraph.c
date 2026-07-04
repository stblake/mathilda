/* icosahedralgraph.c - IcosahedralGraph[]: the graph of the regular icosahedron
 * (a Platonic solid). It is a pentagonal antiprism (two 5-cycle rings joined in
 * antiprism fashion) capped by two apex vertices, one joined to each ring.
 *
 * 12 vertices: top apex = 1, top ring = 2..6, bottom ring = 7..11, bottom apex =
 * 12. Edges (30): the two ring 5-cycles, the apex-to-ring spokes (5 each), and the
 * antiprism links top_r ~ bottom_r and top_r ~ bottom_{r+1}. 5-regular, connected,
 * non-bipartite, chromatic number 4, Hamiltonian. Takes no arguments. O(1).
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_icosahedral_graph(Expr* res) {
    if (res->data.function.arg_count != 0) return NULL;

    enum { V = 12 };
    Expr** vc = malloc(V * sizeof(Expr*));
    for (int i = 0; i < V; i++) vc[i] = expr_new_integer(i + 1);

    Expr** ec = malloc(30 * sizeof(Expr*));
    size_t me = 0;
    #define ADD(a,b) do { Expr* _e[2] = { expr_new_integer(a), expr_new_integer(b) }; \
        ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), _e, 2); } while (0)
    const int TOP = 1, BOT = 12;
    for (int r = 0; r < 5; r++) {
        int t = 2 + r, t1 = 2 + (r + 1) % 5;      /* top ring vertex + next */
        int b = 7 + r, b1 = 7 + (r + 1) % 5;      /* bottom ring vertex + next */
        ADD(t, t1);                                /* top ring cycle */
        ADD(b, b1);                                /* bottom ring cycle */
        ADD(TOP, t);                               /* top apex spoke */
        ADD(BOT, b);                               /* bottom apex spoke */
        ADD(t, b);                                 /* antiprism link */
        ADD(t, b1);                                /* offset antiprism link */
    }
    #undef ADD

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)V);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
