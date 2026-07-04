/* density.c - GraphDensity[g]: the fraction of possible edges that are present,
 * an exact rational in [0, 1].
 *
 *   undirected g -> m / C(n,2) = 2m / (n(n-1))   (n possible unordered pairs)
 *   directed g   -> m / (n(n-1))                 (ordered pairs, no self-loops)
 *
 * A graph with fewer than two vertices has no possible edges; its density is 0
 * by convention. The value is built as num/den and reduced by the evaluator, so
 * it prints as a clean integer or Rational (K_n -> 1, an empty graph -> 0).
 *
 * Memory (SPEC section 4): returns a fresh number; frees res. NULL on a
 * non-graph argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include "eval.h"
#include <stdlib.h>

/* Exact num/den as a reduced value (integer or Rational). den > 0 required. */
static Expr* ratio(long num, long den) {
    Expr* pa[2] = { expr_new_integer(den), expr_new_integer(-1) };
    Expr* inv = expr_new_function(expr_new_symbol(SYM_Power), pa, 2);
    Expr* ta[2] = { expr_new_integer(num), inv };
    return evaluate(expr_new_function(expr_new_symbol(SYM_Times), ta, 2));
}

Expr* builtin_graph_density(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    long n  = (long)verts->data.function.arg_count;
    long ne = (long)edges->data.function.arg_count;
    if (n < 2) return expr_new_integer(0);

    int directed = (ne > 0);
    for (long i = 0; i < ne; i++)
        if (graph_edge_kind(edges->data.function.args[i]) != SYM_DirectedEdge) { directed = 0; break; }

    /* directed: n(n-1) ordered pairs; undirected: n(n-1)/2 unordered pairs. */
    if (directed) return ratio(ne, n * (n - 1));
    return ratio(2 * ne, n * (n - 1));
}
