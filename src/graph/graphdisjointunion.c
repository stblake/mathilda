/* graphdisjointunion.c - GraphDisjointUnion[g1, g2]: the disjoint union of g1 and
 * g2. The vertices are relabelled to 1..n1+n2 (g1's block first, g2's block
 * offset by n1) and the edges are those of g1 and g2 (relabelled, kinds
 * preserved) with NO edges between the two blocks -- so the result has exactly
 * the two graphs as its components.
 *
 * Relabelling to consecutive integers makes the union genuinely disjoint even
 * when the graphs share vertex names. n1+n2 vertices, m1+m2 edges. O(V+E). It is
 * GraphJoin without the n1*n2 cross edges.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL unless
 * both arguments are valid graphs.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Copy an edge with endpoints relabelled to 1-based integers by their index in
 * the source vertex list, offset by `base`. Preserves edge kind. */
static Expr* relabel_edge(const Expr* verts, const Expr* e, int base) {
    const char* kind = graph_edge_kind(e);
    int ia = graph_vertex_index(verts, e->data.function.args[0]);
    int ib = graph_vertex_index(verts, e->data.function.args[1]);
    Expr* args[2] = { expr_new_integer(base + ia + 1), expr_new_integer(base + ib + 1) };
    return expr_new_function(expr_new_symbol(kind), args, 2);
}

Expr* builtin_graph_disjoint_union(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g1 = res->data.function.args[0];
    const Expr* g2 = res->data.function.args[1];
    if (!graph_is_valid(g1) || !graph_is_valid(g2)) return NULL;

    const Expr* v1 = g1->data.function.args[0];
    const Expr* e1 = g1->data.function.args[1];
    const Expr* v2 = g2->data.function.args[0];
    const Expr* e2 = g2->data.function.args[1];
    int n1 = (int)v1->data.function.arg_count;
    int n2 = (int)v2->data.function.arg_count;
    size_t m1 = e1->data.function.arg_count, m2 = e2->data.function.arg_count;

    int N = n1 + n2;
    Expr** vc = malloc((N > 0 ? N : 1) * sizeof(Expr*));
    for (int i = 0; i < N; i++) vc[i] = expr_new_integer(i + 1);

    size_t total = m1 + m2;
    Expr** ec = (total > 0) ? malloc(total * sizeof(Expr*)) : NULL;
    size_t me = 0;
    for (size_t i = 0; i < m1; i++) ec[me++] = relabel_edge(v1, e1->data.function.args[i], 0);
    for (size_t i = 0; i < m2; i++) ec[me++] = relabel_edge(v2, e2->data.function.args[i], n1);

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)N);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
