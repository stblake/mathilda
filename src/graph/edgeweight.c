/* edgeweight.c - EdgeWeight[g]: the graph's per-edge weights, in EdgeList
 * order. Defaults to List[1, 1, ..., 1] (one per edge) when g carries no
 * EdgeWeight -- matching Wolfram Language's own behavior for an unweighted
 * graph, and giving WeightedAdjacencyMatrix[g] a well-defined answer for
 * every valid graph, not just ones explicitly built with weights.
 *
 * Memory (SPEC section 4): returns a fresh list; the evaluator frees res.
 */

#include "graph.h"
#include "expr.h"

Expr* builtin_edge_weight(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    return graph_resolve_edge_weights(g);   /* NULL (unevaluated) if g is not a valid graph */
}
