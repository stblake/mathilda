/* turangraph.c - TuranGraph[n, r]: the Turan graph T(n, r), the complete
 * r-partite graph on n vertices whose parts are as equal in size as possible.
 * Two vertices are adjacent iff they lie in different parts, so it is the
 * n-vertex graph with the most edges that contains no (r+1)-clique (Turan's
 * theorem).
 *
 * Vertices are 1..n; vertex i (0-based) is placed in part i mod r, which yields
 * balanced parts (sizes differ by at most 1). Undirected edges join every pair
 * in different parts. O(n^2). T(n, 1) is edgeless, T(n, n) is K_n, T(4, 2) = C4,
 * T(6, 3) is the octahedron.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL unless
 * n >= 0 and r >= 1 are integers.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_turan_graph(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* ne = res->data.function.args[0];
    const Expr* re = res->data.function.args[1];
    if (ne->type != EXPR_INTEGER || re->type != EXPR_INTEGER) return NULL;
    long n = (long)ne->data.integer, r = (long)re->data.integer;
    if (n < 0 || r < 1) return NULL;

    Expr** vc = malloc((n > 0 ? (size_t)n : 1) * sizeof(Expr*));
    for (long i = 0; i < n; i++) vc[i] = expr_new_integer(i + 1);

    /* Upper bound on edges: n(n-1)/2. */
    size_t cap = (size_t)n * (n > 0 ? (size_t)(n - 1) : 0) / 2;
    Expr** ec = (cap > 0) ? malloc(cap * sizeof(Expr*)) : NULL;
    size_t me = 0;
    for (long i = 0; i < n; i++)
        for (long j = i + 1; j < n; j++)
            if (i % r != j % r) {                  /* different parts */
                Expr* args[2] = { expr_new_integer(i + 1), expr_new_integer(j + 1) };
                ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), args, 2);
            }

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)n);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
