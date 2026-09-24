/* membership.c - VertexQ[g, v] and EdgeQ[g, e].
 *
 *   VertexQ[g, v]   True iff v is a vertex of g
 *   EdgeQ[g, e]     True iff e is an edge of g
 *
 * Membership is structural (SameQ, via expr_eq), as in the Wolfram Language:
 * a vertex 1 is not matched by 1.0.
 *
 * EdgeQ accepts the same edge sugar the constructor does -- u -> v means
 * DirectedEdge[u, v], u <-> v means UndirectedEdge[u, v]. An undirected query
 * matches an UndirectedEdge in either orientation; a directed query matches
 * only a DirectedEdge with the same ordered endpoints. Direction is never
 * blurred: u -> v is not an edge of Graph[{u, v}, {u <-> v}].
 *
 * Both give False for anything that is not a valid graph. Both are O(1) hash
 * probes into the validated-graph memo (graph_util.c), which already holds the
 * vertex index and edge-key set that validating g built.
 *
 * Memory (SPEC section 4): returns a fresh symbol; the evaluator frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"

static Expr* truth(int b) { return expr_new_symbol(b ? SYM_True : SYM_False); }

Expr* builtin_vertex_q(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    return truth(graph_vertex_position(res->data.function.args[0],
                                       res->data.function.args[1]) >= 0);
}

/* The canonical edge head a query denotes (DirectedEdge/UndirectedEdge), or
 * NULL when q is not a 2-argument edge or edge sugar. */
static const char* query_kind(const Expr* q) {
    if (!q || q->type != EXPR_FUNCTION || q->data.function.arg_count != 2) return NULL;
    const Expr* h = q->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return NULL;
    const char* name = h->data.symbol.name;
    if (name == SYM_DirectedEdge || name == SYM_Rule) return SYM_DirectedEdge;
    if (name == SYM_UndirectedEdge || name == SYM_TwoWayRule) return SYM_UndirectedEdge;
    return NULL;
}

Expr* builtin_edge_q(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return truth(0);
    const Expr* q = res->data.function.args[1];
    const char* kind = query_kind(q);
    if (!kind) return truth(0);
    return truth(graph_has_edge(g, q->data.function.args[0], q->data.function.args[1],
                                kind == SYM_DirectedEdge) > 0);
}
