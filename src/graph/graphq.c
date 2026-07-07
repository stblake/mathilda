/* graphq.c - GraphQ[g]: is g a valid graph?
 *
 * A thin wrapper over graph_is_valid (graph_util.c): returns the symbol True
 * when the (already-evaluated) argument is a canonical, valid graph, and False
 * otherwise. Non-unary calls are left unevaluated (NULL).
 *
 * Memory (SPEC section 4): returns a freshly-allocated symbol; the evaluator
 * frees `res`.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"

Expr* builtin_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* arg = res->data.function.args[0];
    return expr_new_symbol(graph_is_valid(arg) ? SYM_True : SYM_False);
}
