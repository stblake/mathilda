/* completekarytree.c - CompleteKaryTree[L] / CompleteKaryTree[L, k]: the complete
 * k-ary tree with L levels (k = 2 by default) -- a root whose every internal node
 * has exactly k children, filled to depth L.
 *
 * Vertices are 1..V in breadth-first (heap) order, so vertex i >= 2 has parent
 * floor((i-2)/k) + 1; the V-1 tree edges are (parent(i), i), undirected. The
 * vertex count is V = (k^L - 1)/(k - 1) for k >= 2 and V = L for k = 1 (a path).
 * O(V). The result is always a tree (TreeGraphQ is True).
 *
 * L >= 1 and k >= 1 must be integers; the call is left unevaluated if the tree
 * would be larger than a safety cap.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL on bad
 * arguments.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

#define KARY_MAX_VERTICES 1000000L

Expr* builtin_complete_kary_tree(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    const Expr* Le = res->data.function.args[0];
    if (Le->type != EXPR_INTEGER) return NULL;
    long L = (long)Le->data.integer;
    long k = 2;
    if (argc == 2) {
        const Expr* ke = res->data.function.args[1];
        if (ke->type != EXPR_INTEGER) return NULL;
        k = (long)ke->data.integer;
    }
    if (L < 1 || k < 1) return NULL;

    /* Vertex count V, with overflow/size guard. */
    long V;
    if (k == 1) {
        V = L;
    } else {
        V = 1; long level = 1;
        for (long d = 1; d < L; d++) {
            level *= k;                                 /* nodes at depth d */
            if (level > KARY_MAX_VERTICES || V > KARY_MAX_VERTICES - level) return NULL;
            V += level;
        }
    }
    if (V > KARY_MAX_VERTICES) return NULL;

    Expr** vc = malloc((size_t)V * sizeof(Expr*));
    for (long i = 0; i < V; i++) vc[i] = expr_new_integer(i + 1);

    Expr** ec = (V > 1) ? malloc((size_t)(V - 1) * sizeof(Expr*)) : NULL;
    size_t me = 0;
    for (long i = 2; i <= V; i++) {
        long parent = (i - 2) / k + 1;
        Expr* args[2] = { expr_new_integer(parent), expr_new_integer(i) };
        ec[me++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), args, 2);
    }

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)V);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
