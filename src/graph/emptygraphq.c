/* emptygraphq.c - EmptyGraphQ[g]: True iff g has no edges (an edgeless graph on
 * any number of vertices, including none). O(1) -- just the edge count.
 *
 * Memory (SPEC section 4): returns True/False; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"

Expr* builtin_empty_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    int empty = (g->data.function.args[1]->data.function.arg_count == 0);
    return expr_new_symbol(empty ? SYM_True : SYM_False);
}
