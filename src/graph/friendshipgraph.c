/* friendshipgraph.c - FriendshipGraph[n]: the friendship (windmill) graph F_n,
 * n triangles all sharing a single common vertex.
 *
 * Vertex 1 is the hub; the outer vertices 2..2n+1 are paired (2+2i, 3+2i), and
 * each pair forms a triangle with the hub. So there are 2n+1 vertices and 3n
 * edges: for each pair, the two spokes hub-a, hub-b and the rim edge a-b. O(n).
 * F_1 is the triangle K_3, F_2 the bowtie. By the friendship theorem it is the
 * unique graph in which every two vertices have exactly one common neighbour.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL unless
 * n >= 1 is an integer.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_friendship_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* ne = res->data.function.args[0];
    if (ne->type != EXPR_INTEGER) return NULL;
    long n = (long)ne->data.integer;
    if (n < 1) return NULL;

    long V = 2 * n + 1;
    Expr** vc = malloc((size_t)V * sizeof(Expr*));
    for (long i = 0; i < V; i++) vc[i] = expr_new_integer(i + 1);

    Expr** ec = malloc((size_t)(3 * n) * sizeof(Expr*));
    size_t me = 0;
    for (long i = 0; i < n; i++) {
        long a = 2 + 2 * i, b = 3 + 2 * i;         /* a triangle's rim pair */
        long tri[3][2] = { {1, a}, {1, b}, {a, b} };
        for (int e = 0; e < 3; e++) {
            Expr* args[2] = { expr_new_integer(tri[e][0]), expr_new_integer(tri[e][1]) };
            ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), args, 2);
        }
    }

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)V);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
