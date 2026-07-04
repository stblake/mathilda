/* cocktailparty.c - CocktailPartyGraph[n]: the cocktail-party graph K_{n x 2},
 * the complete n-partite graph with all n parts of size 2. Equivalently, 2n
 * guests in n couples where everyone talks to everyone except their own partner
 * -- the complement of a perfect matching on 2n vertices.
 *
 * Vertices are 1..2n; vertex i (0-based) is in couple i/2, and two vertices are
 * joined iff they belong to different couples. It is (2n-2)-regular with
 * 2n(n-1) edges. O(n^2). CocktailPartyGraph[2] = C4, CocktailPartyGraph[3] is the
 * octahedron.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL unless
 * n >= 1 is an integer.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_cocktail_party_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* ne = res->data.function.args[0];
    if (ne->type != EXPR_INTEGER) return NULL;
    long n = (long)ne->data.integer;
    if (n < 1) return NULL;

    long V = 2 * n;
    Expr** vc = malloc((size_t)V * sizeof(Expr*));
    for (long i = 0; i < V; i++) vc[i] = expr_new_integer(i + 1);

    size_t cap = (size_t)(2 * n * (n - 1));            /* 2n(n-1) edges */
    Expr** ec = (cap > 0) ? malloc(cap * sizeof(Expr*)) : NULL;
    size_t me = 0;
    for (long i = 0; i < V; i++)
        for (long j = i + 1; j < V; j++)
            if (i / 2 != j / 2) {                       /* different couples */
                Expr* a[2] = { expr_new_integer(i + 1), expr_new_integer(j + 1) };
                ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), a, 2);
            }

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)V);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
