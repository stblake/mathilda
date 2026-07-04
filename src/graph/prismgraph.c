/* prismgraph.c - PrismGraph[n]: the n-gonal prism (a.k.a. circular ladder), two
 * concentric n-cycles joined rung-by-rung -- the Cartesian product C_n [] K_2.
 *
 * Vertices are the outer cycle 1..n and the inner cycle n+1..2n. Edges are the
 * outer cycle o_i ~ o_{i+1}, the inner cycle c_i ~ c_{i+1}, and the rungs
 * o_i ~ c_i (indices mod n). 3-regular with 3n edges. O(n). PrismGraph[3] is the
 * triangular prism, PrismGraph[4] the cube. (Isomorphic to
 * GeneralizedPetersenGraph[n, 1], provided under its own common name.)
 *
 * n >= 3 must be an integer. Memory (SPEC section 4): returns a freshly-built
 * Graph; frees res. NULL otherwise.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_prism_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* ne = res->data.function.args[0];
    if (ne->type != EXPR_INTEGER) return NULL;
    long n = (long)ne->data.integer;
    if (n < 3) return NULL;

    long V = 2 * n;
    Expr** vc = malloc((size_t)V * sizeof(Expr*));
    for (long i = 0; i < V; i++) vc[i] = expr_new_integer(i + 1);

    Expr** ec = malloc((size_t)(3 * n) * sizeof(Expr*));
    size_t me = 0;
    #define ADD(a,b) do { Expr* _e[2] = { expr_new_integer(a), expr_new_integer(b) }; \
        ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), _e, 2); } while (0)
    for (long i = 0; i < n; i++) {
        long o = i + 1, o1 = (i + 1) % n + 1, c = n + i + 1, c1 = n + (i + 1) % n + 1;
        ADD(o, o1);          /* outer cycle */
        ADD(c, c1);          /* inner cycle */
        ADD(o, c);           /* rung */
    }
    #undef ADD

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)V);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
