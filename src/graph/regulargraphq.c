/* regulargraphq.c - RegularGraphQ[g]: True iff g is regular -- every vertex has
 * the same degree. For a directed graph this means all out-degrees are equal
 * and all in-degrees are equal (an undirected graph has out == in per vertex,
 * so it reduces to "all degrees equal").
 *
 * One pass over the adjacency degree arrays, O(V). A graph with 0 or 1 vertex is
 * vacuously regular.
 *
 * Memory (SPEC section 4): returns True/False; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_regular_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GraphAdj* a = graph_build_adj(res->data.function.args[0]);
    if (!a) return NULL;
    int n = a->n;

    int regular = 1;
    for (int i = 1; i < n; i++)
        if (a->outdeg[i] != a->outdeg[0] || a->indeg[i] != a->indeg[0]) { regular = 0; break; }
    graph_adj_free(a);
    return expr_new_symbol(regular ? SYM_True : SYM_False);
}
