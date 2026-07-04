/* vertexdelete.c - VertexDelete[g, v] / VertexDelete[g, {v1, ...}]: g with the
 * given vertex (or vertices) removed, along with every edge incident to a removed
 * vertex.
 *
 * A List second argument names several vertices to delete (matching Wolfram's
 * convention); any other expression is a single vertex. Surviving vertices keep
 * their order and edge kinds are preserved. O(V + E) with linear membership
 * tests. VertexDelete[K4, 1] is K3; VertexDelete[PathGraph[3], 2] leaves two
 * isolated vertices.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL on a
 * non-graph first argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Is vertex-expression v named by the delete spec? */
static int in_delete_spec(const Expr* spec, const Expr* v) {
    if (graph_is_list(spec)) {
        for (size_t i = 0; i < spec->data.function.arg_count; i++)
            if (expr_eq(spec->data.function.args[i], v)) return 1;
        return 0;
    }
    return expr_eq(spec, v);
}

Expr* builtin_vertex_delete(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    const Expr* spec = res->data.function.args[1];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t m = edges->data.function.arg_count;

    char* removed = (n > 0) ? calloc((size_t)n, 1) : NULL;
    Expr** vc = malloc((size_t)(n > 0 ? n : 1) * sizeof(Expr*));
    size_t nv = 0;
    for (int i = 0; i < n; i++) {
        const Expr* v = verts->data.function.args[i];
        if (in_delete_spec(spec, v)) { removed[i] = 1; continue; }
        vc[nv++] = expr_copy((Expr*)v);
    }

    Expr** ec = (m > 0) ? malloc(m * sizeof(Expr*)) : NULL;
    size_t me = 0;
    for (size_t k = 0; k < m; k++) {
        const Expr* e = edges->data.function.args[k];
        int a = graph_vertex_index(verts, e->data.function.args[0]);
        int b = graph_vertex_index(verts, e->data.function.args[1]);
        if (a >= 0 && removed[a]) continue;
        if (b >= 0 && removed[b]) continue;
        ec[me++] = expr_copy((Expr*)e);
    }
    free(removed);

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, nv);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
