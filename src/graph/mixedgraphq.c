/* mixedgraphq.c - MixedGraphQ[g]: True iff g has at least one directed edge and
 * at least one undirected edge (a mixed graph). A purely directed, purely
 * undirected, or edgeless graph is not mixed. O(E).
 *
 * Memory (SPEC section 4): returns True/False; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"

Expr* builtin_mixed_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* edges = g->data.function.args[1];
    int has_dir = 0, has_undir = 0;
    for (size_t i = 0; i < edges->data.function.arg_count; i++) {
        const char* k = graph_edge_kind(edges->data.function.args[i]);
        if (k == SYM_DirectedEdge) has_dir = 1;
        else if (k == SYM_UndirectedEdge) has_undir = 1;
        if (has_dir && has_undir) break;
    }
    return expr_new_symbol((has_dir && has_undir) ? SYM_True : SYM_False);
}
