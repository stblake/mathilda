/* edgeadd.c - EdgeAdd[g, e] / EdgeAdd[g, {e1, ...}]: g with the given edge(s)
 * added. Any endpoint that is not already a vertex of g is added as a new vertex,
 * so the result is always a valid graph.
 *
 * An edge spec may be DirectedEdge/UndirectedEdge or the sugar a->b (Rule) /
 * a<->b (TwoWayRule); it is normalized to the corresponding edge kind. To keep
 * the simple-graph invariant, a self-loop (a == b) is skipped and an edge equal
 * to one already present (edge-kind-aware, symmetric for undirected) is not
 * duplicated. O((V+E) * #specs) at small-graph scale.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL on a
 * non-graph first argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

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

static int same_edge(const Expr* ge, const char* kind, const Expr* sa, const Expr* sb) {
    if (graph_edge_kind(ge) != kind) return 0;
    const Expr* a = ge->data.function.args[0];
    const Expr* b = ge->data.function.args[1];
    if (expr_eq(a, sa) && expr_eq(b, sb)) return 1;
    if (kind == SYM_UndirectedEdge && expr_eq(a, sb) && expr_eq(b, sa)) return 1;
    return 0;
}

Expr* builtin_edge_add(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    const Expr* spec = res->data.function.args[1];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    size_t n = verts->data.function.arg_count;
    size_t m = edges->data.function.arg_count;

    /* Spec edges: either one edge or a List of them. */
    const Expr* const* specs; size_t ns; const Expr* one[1];
    if (graph_is_list(spec)) { specs = (const Expr* const*)spec->data.function.args; ns = spec->data.function.arg_count; }
    else { one[0] = spec; specs = one; ns = 1; }

    /* Vertex pointers (borrowed): g's vertices plus any new endpoints. */
    const Expr** vptr = malloc((n + 2 * ns) * sizeof(Expr*));
    size_t nv = 0;
    for (size_t i = 0; i < n; i++) vptr[nv++] = verts->data.function.args[i];

    /* Output edges (owned copies): g's edges plus accepted new ones. */
    Expr** ec = malloc((m + ns) * sizeof(Expr*));
    size_t me = 0;
    for (size_t i = 0; i < m; i++) ec[me++] = expr_copy(edges->data.function.args[i]);

    for (size_t s = 0; s < ns; s++) {
        const char* kind; const Expr *a, *b;
        if (!spec_edge(specs[s], &kind, &a, &b)) continue;
        if (expr_eq(a, b)) continue;                       /* no self-loops */
        int dup = 0;
        for (size_t j = 0; j < me; j++) if (same_edge(ec[j], kind, a, b)) { dup = 1; break; }
        if (dup) continue;
        for (int which = 0; which < 2; which++) {          /* register endpoints */
            const Expr* v = which ? b : a;
            int seen = 0;
            for (size_t j = 0; j < nv; j++) if (expr_eq(vptr[j], v)) { seen = 1; break; }
            if (!seen) vptr[nv++] = v;
        }
        Expr* args[2] = { expr_copy((Expr*)a), expr_copy((Expr*)b) };
        ec[me++] = expr_new_function(expr_new_symbol(kind), args, 2);
    }

    Expr** vc = malloc((nv > 0 ? nv : 1) * sizeof(Expr*));
    for (size_t i = 0; i < nv; i++) vc[i] = expr_copy((Expr*)vptr[i]);
    free(vptr);

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, nv);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
