/* reciprocity.c - GraphReciprocity[g]: the fraction of arcs whose reverse is
 * also present, an exact rational in [0, 1]. It measures how mutual a directed
 * graph's connections are; an undirected graph is fully reciprocal (1).
 *
 * Model everything as a set of directed arcs on an n x n boolean matrix: a
 * DirectedEdge a->b sets one arc, an UndirectedEdge {a,b} sets both a->b and
 * b->a (an undirected edge is inherently mutual). Then
 *
 *   reciprocity = #{arcs (a,b) with (b,a) also present} / #arcs.
 *
 * For a purely directed graph this is the usual "fraction of edges that are
 * reciprocated"; for an undirected graph every arc has its reverse, giving 1.
 * With no edges the ratio is undefined, reported as 0 by convention. O(V^2).
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

Expr* builtin_graph_reciprocity(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    char* arc = (n > 0) ? calloc((size_t)n * (size_t)n, 1) : NULL;
    for (size_t k = 0; k < ne; k++) {
        const Expr* ed = edges->data.function.args[k];
        int ia = graph_vertex_index(verts, ed->data.function.args[0]);
        int ib = graph_vertex_index(verts, ed->data.function.args[1]);
        if (ia < 0 || ib < 0 || ia == ib) continue;
        arc[(size_t)ia * n + ib] = 1;
        if (graph_edge_kind(ed) == SYM_UndirectedEdge) arc[(size_t)ib * n + ia] = 1;
    }

    long total = 0, recip = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (arc[(size_t)i * n + j]) {
                total++;
                if (arc[(size_t)j * n + i]) recip++;
            }
    free(arc);

    if (total == 0) return expr_new_integer(0);   /* no arcs -> undefined, 0 */
    return ratio(recip, total);
}
