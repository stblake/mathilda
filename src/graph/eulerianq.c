/* eulerianq.c - EulerianGraphQ[g]: True iff g has an Eulerian cycle (a closed
 * walk traversing every edge exactly once).
 *
 * Conditions:
 *   undirected: every vertex has even degree, and all vertices of nonzero degree
 *               lie in one connected component;
 *   directed:   in-degree equals out-degree at every vertex, and all nonzero-
 *               degree vertices lie in one component. (When in==out everywhere,
 *               weak connectivity implies strong connectivity, so the undirected
 *               component test suffices.)
 * Isolated vertices are ignored; an edgeless graph is vacuously Eulerian.
 *
 * O(V+E) — degree scan plus one component count over the shared adjacency.
 * Memory (SPEC section 4): returns a fresh True/False symbol; frees res. False
 * on a non-graph argument (matching GraphQ).
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_eulerian_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return expr_new_symbol(SYM_False);

    const Expr* edges = g->data.function.args[1];
    size_t ne = edges->data.function.arg_count;
    int directed = (ne > 0);
    for (size_t i = 0; i < ne; i++)
        if (graph_edge_kind(edges->data.function.args[i]) != SYM_DirectedEdge) { directed = 0; break; }

    GraphAdj* a = graph_build_adj(g);
    if (!a) return expr_new_symbol(SYM_False);
    int n = a->n;

    /* Degree condition. For an undirected graph a->outdeg is the (incident)
     * degree; for a directed graph we need in-degree == out-degree. */
    int degree_ok = 1;
    for (int i = 0; i < n; i++) {
        if (directed) { if (a->indeg[i] != a->outdeg[i]) { degree_ok = 0; break; } }
        else          { if (a->outdeg[i] % 2 != 0)       { degree_ok = 0; break; } }
    }

    int eulerian = 0;
    if (degree_ok) {
        /* Connectivity over the nonzero-degree vertices only. */
        char* removed = (n > 0) ? calloc((size_t)n, 1) : NULL;
        for (int i = 0; i < n; i++)
            if (a->outdeg[i] + a->indeg[i] == 0) removed[i] = 1;   /* isolated */
        int active = 0;
        int comps = graph_count_components(a, removed, &active);
        eulerian = (active == 0) ? 1 : (comps == 1);   /* edgeless is vacuous */
        free(removed);
    }

    graph_adj_free(a);
    return expr_new_symbol(eulerian ? SYM_True : SYM_False);
}
