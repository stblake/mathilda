/* antiprismgraph.c - AntiprismGraph[n]: the n-antiprism, two concentric n-cycles
 * (outer 1..n, inner n+1..2n) whose rings are offset so each outer vertex joins
 * two consecutive inner vertices.
 *
 * Edges: the outer cycle o_i ~ o_{i+1}, the inner cycle c_i ~ c_{i+1}, and the
 * cross edges o_i ~ c_i and o_i ~ c_{i+1} (indices mod n). Every vertex has
 * degree 4 (two ring + two cross), giving a 4-regular graph with 4n edges. O(n^2)
 * via a boolean adjacency. AntiprismGraph[3] is the octahedron.
 *
 * n >= 3 must be an integer. Memory (SPEC section 4): returns a freshly-built
 * Graph; frees res. NULL otherwise.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_antiprism_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* ne = res->data.function.args[0];
    if (ne->type != EXPR_INTEGER) return NULL;
    long n = (long)ne->data.integer;
    if (n < 3) return NULL;

    long V = 2 * n;
    char* adj = calloc((size_t)V * (size_t)V, 1);
    #define SET(a,b) do { adj[(size_t)(a) * V + (b)] = adj[(size_t)(b) * V + (a)] = 1; } while (0)
    for (long i = 0; i < n; i++) {
        long o = i, c = n + i, o1 = (i + 1) % n, c1 = n + (i + 1) % n;
        SET(o, o1);          /* outer cycle */
        SET(c, c1);          /* inner cycle */
        SET(o, c);           /* cross */
        SET(o, c1);          /* offset cross */
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
