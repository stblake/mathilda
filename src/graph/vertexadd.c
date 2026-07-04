/* vertexadd.c - VertexAdd[g, v] / VertexAdd[g, {v1, ...}]: g with the given
 * vertex (or vertices) added as isolated vertices. The edge set is unchanged.
 *
 * A List second argument names several vertices to add (Wolfram's convention);
 * any other expression is a single vertex. A vertex already present is not
 * duplicated. O((V + #new) * #new). New vertices are appended after the existing
 * ones.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL on a
 * non-graph first argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_vertex_add(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    const Expr* spec = res->data.function.args[1];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    size_t n = verts->data.function.arg_count;
    size_t m = edges->data.function.arg_count;

    /* New vertices: a List names several, anything else is a single vertex. */
    const Expr* const* adds; size_t na; const Expr* one[1];
    if (graph_is_list(spec)) { adds = (const Expr* const*)spec->data.function.args; na = spec->data.function.arg_count; }
    else { one[0] = spec; adds = one; na = 1; }

    Expr** vc = malloc((n + na) * sizeof(Expr*));
    size_t nv = 0;
    for (size_t i = 0; i < n; i++) vc[nv++] = expr_copy(verts->data.function.args[i]);
    for (size_t i = 0; i < na; i++) {
        const Expr* v = adds[i];
        int seen = 0;
        for (size_t j = 0; j < nv; j++) if (expr_eq(vc[j], v)) { seen = 1; break; }
        if (!seen) vc[nv++] = expr_copy((Expr*)v);
    }

    Expr** ec = (m > 0) ? malloc(m * sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < m; i++) ec[i] = expr_copy(edges->data.function.args[i]);

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, nv);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, m);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
