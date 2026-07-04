/* dodecahedralgraph.c - DodecahedralGraph[]: the graph of the regular
 * dodecahedron (one of the five Platonic solids), which is exactly the
 * generalized Petersen graph GP(10, 2).
 *
 * 20 vertices (outer 10-cycle 1..10, inner set 11..20) and 30 edges: the outer
 * cycle o_j ~ o_{j+1}, the spokes o_j ~ i_j, and the inner star polygon
 * i_j ~ i_{j+2} (indices mod 10). 3-regular, girth 5, non-bipartite, connected.
 * Takes no arguments. O(1).
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_dodecahedral_graph(Expr* res) {
    if (res->data.function.arg_count != 0) return NULL;

    enum { N = 10, V = 20 };
    Expr** vc = malloc(V * sizeof(Expr*));
    for (int i = 0; i < V; i++) vc[i] = expr_new_integer(i + 1);

    Expr** ec = malloc(30 * sizeof(Expr*));
    size_t me = 0;
    #define ADD(a,b) do { Expr* _e[2] = { expr_new_integer(a), expr_new_integer(b) }; \
        ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), _e, 2); } while (0)
    for (int j = 0; j < N; j++) {
        ADD(j + 1, (j + 1) % N + 1);            /* outer cycle */
        ADD(j + 1, N + 1 + j);                  /* spoke */
        ADD(N + 1 + j, N + 1 + (j + 2) % N);    /* inner star polygon {10/2} */
    }
    #undef ADD

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)V);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
