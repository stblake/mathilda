/* generalizedpetersen.c - GeneralizedPetersenGraph[n, k]: the generalized
 * Petersen graph GP(n, k). It has 2n vertices -- an outer n-cycle 1..n and an
 * inner set n+1..2n -- with three kinds of edges: the outer cycle o_j ~ o_{j+1},
 * the spokes o_j ~ i_j, and the inner "star polygon" i_j ~ i_{j+k} (indices mod
 * n).
 *
 * A 2n x 2n boolean adjacency dedups the inner edges (when 2k = n they form a
 * matching rather than a cycle). O(n^2). Famous members: GP(n,1) is the n-prism
 * (GP(4,1) the cube), GP(5,2) is the Petersen graph, GP(6,2) the Durer graph,
 * GP(8,3) the Mobius-Kantor graph, GP(10,3) the Desargues graph.
 *
 * n >= 3 and 1 <= k < n must be integers.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL on bad
 * arguments.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_generalized_petersen_graph(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* ne = res->data.function.args[0];
    const Expr* ke = res->data.function.args[1];
    if (ne->type != EXPR_INTEGER || ke->type != EXPR_INTEGER) return NULL;
    long n = (long)ne->data.integer, k = (long)ke->data.integer;
    if (n < 3 || k < 1 || k >= n) return NULL;

    long V = 2 * n;
    char* adj = calloc((size_t)V * (size_t)V, 1);
    #define SET(a,b) do { adj[(size_t)(a) * V + (b)] = adj[(size_t)(b) * V + (a)] = 1; } while (0)
    for (long j = 0; j < n; j++) {
        long o = j, oi = n + j;
        SET(o, (j + 1) % n);                 /* outer cycle */
        SET(o, oi);                          /* spoke */
        long inner = n + (j + k) % n;        /* inner star polygon */
        if (inner != oi) SET(oi, inner);
    }
    #undef SET

    Expr** vc = malloc((size_t)V * sizeof(Expr*));
    for (long i = 0; i < V; i++) vc[i] = expr_new_integer(i + 1);

    size_t cap = (size_t)V * (size_t)(V - 1) / 2;
    Expr** ec = malloc(cap * sizeof(Expr*));
    size_t me = 0;
    for (long i = 0; i < V; i++)
        for (long j = i + 1; j < V; j++)
            if (adj[(size_t)i * V + j]) {
                Expr* a[2] = { expr_new_integer(i + 1), expr_new_integer(j + 1) };
                ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), a, 2);
            }
    free(adj);

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)V);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
