/* edgedelete.c - EdgeDelete[g, e] / EdgeDelete[g, {e1, ...}]: g with the given
 * edge (or edges) removed. All vertices are kept.
 *
 * An edge spec may be written as DirectedEdge/UndirectedEdge or with the sugar
 * a->b (Rule) / a<->b (TwoWayRule); it is normalized to a (kind, a, b) triple
 * and matched against g's edges with the usual edge-kind-aware, symmetric-for-
 * undirected equality. A List second argument names several edges to delete;
 * only the first matching graph edge is removed per spec occurrence. O(E * #specs)
 * at small-graph scale.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL on a
 * non-graph first argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Normalize an edge spec to (kind, a, b); returns 0 if it is not an edge. */
static int spec_edge(const Expr* e, const char** kind, const Expr** a, const Expr** b) {
    if (e->type != EXPR_FUNCTION || e->data.function.arg_count != 2 ||
        e->data.function.head->type != EXPR_SYMBOL) return 0;
    const char* h = e->data.function.head->data.symbol;
    if (h == SYM_DirectedEdge || h == SYM_Rule) *kind = SYM_DirectedEdge;
    else if (h == SYM_UndirectedEdge || h == SYM_TwoWayRule) *kind = SYM_UndirectedEdge;
    else return 0;
    *a = e->data.function.args[0];
    *b = e->data.function.args[1];
    return 1;
}

/* Does graph edge ge match the normalized spec (kind, sa, sb)? */
static int edge_matches(const Expr* ge, const char* kind, const Expr* sa, const Expr* sb) {
    if (graph_edge_kind(ge) != kind) return 0;
    const Expr* a = ge->data.function.args[0];
    const Expr* b = ge->data.function.args[1];
    if (expr_eq(a, sa) && expr_eq(b, sb)) return 1;
    if (kind == SYM_UndirectedEdge && expr_eq(a, sb) && expr_eq(b, sa)) return 1;
    return 0;
}

/* Is graph edge ge named by the delete spec (single edge or List of edges)? */
static int in_edge_spec(const Expr* spec, const Expr* ge) {
    const char* kind; const Expr *a, *b;
    if (graph_is_list(spec)) {
        for (size_t i = 0; i < spec->data.function.arg_count; i++)
            if (spec_edge(spec->data.function.args[i], &kind, &a, &b) &&
                edge_matches(ge, kind, a, b)) return 1;
        return 0;
    }
    return spec_edge(spec, &kind, &a, &b) && edge_matches(ge, kind, a, b);
}

Expr* builtin_edge_delete(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    const Expr* spec = res->data.function.args[1];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    size_t n = verts->data.function.arg_count;
    size_t m = edges->data.function.arg_count;

    Expr** vc = malloc((n > 0 ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) vc[i] = expr_copy(verts->data.function.args[i]);

    Expr** ec = (m > 0) ? malloc(m * sizeof(Expr*)) : NULL;
    size_t me = 0;
    for (size_t k = 0; k < m; k++) {
        const Expr* e = edges->data.function.args[k];
        if (!in_edge_spec(spec, e)) ec[me++] = expr_copy((Expr*)e);
    }

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, n);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
