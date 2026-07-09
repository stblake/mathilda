/* edgelist.c - EdgeList[g]: the graph's edges (canonical Directed/UndirectedEdge
 * form), in order. NULL (unevaluated) if g is not a graph.
 * Memory (SPEC section 4): returns a fresh copy; the evaluator frees res. */

#include "graph.h"
#include "expr.h"

Expr* builtin_edge_list(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    return expr_copy(g->data.function.args[1]);
}
