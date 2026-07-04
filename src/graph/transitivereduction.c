/* transitivereduction.c - TransitiveReductionGraph[g]: the transitive reduction
 * of a directed acyclic graph -- the graph on the same vertices with the fewest
 * edges that preserves reachability. For a DAG this is unique: keep edge u->v
 * iff there is no length->=2 path from u to v, i.e. no intermediate w with
 * u ~> w ~> v.
 *
 * Reachability (paths of length >= 1) is found by a BFS from each vertex over
 * the successor adjacency; self-reachability signals a cycle. Transitive
 * reduction is only well defined (and preserving) on acyclic digraphs, so the
 * call is left unevaluated when g has a directed cycle (an undirected edge is a
 * 2-cycle, so any undirected edge also declines). O(V*(V+E) + E*V).
 *
 * Memory (SPEC section 4): returns a freshly-built Graph; frees res. NULL on a
 * non-graph or a graph with a cycle.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_transitive_reduction_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    int n = a->n;

    /* reach[u*n+v] = 1 iff v is reachable from u by a path of length >= 1. */
    char* reach = (n > 0) ? calloc((size_t)n * (size_t)n, 1) : NULL;
    int* q = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;
    int cyclic = 0;
    for (int s = 0; s < n; s++) {
        char* row = reach + (size_t)s * n;
        int head = 0, tail = 0;
        for (int j = 0; j < a->outdeg[s]; j++) {
            int w = a->out[s][j];
            if (!row[w]) { row[w] = 1; q[tail++] = w; }
        }
        while (head < tail) {
            int v = q[head++];
            for (int j = 0; j < a->outdeg[v]; j++) {
                int w = a->out[v][j];
                if (!row[w]) { row[w] = 1; q[tail++] = w; }
            }
        }
        if (row[s]) { cyclic = 1; break; }        /* s reaches itself -> cycle */
    }
    free(q);
    if (cyclic) { free(reach); graph_adj_free(a); return NULL; }

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    size_t m = edges->data.function.arg_count;

    Expr** vc = malloc((size_t)(n > 0 ? n : 1) * sizeof(Expr*));
    for (int i = 0; i < n; i++) vc[i] = expr_copy(verts->data.function.args[i]);

    Expr** ec = (m > 0) ? malloc(m * sizeof(Expr*)) : NULL;
    size_t me = 0;
    for (size_t k = 0; k < m; k++) {
        const Expr* e = edges->data.function.args[k];
        int u = graph_vertex_index(verts, e->data.function.args[0]);
        int v = graph_vertex_index(verts, e->data.function.args[1]);
        int redundant = 0;
        for (int w = 0; w < n && !redundant; w++)
            if (w != u && w != v && reach[(size_t)u * n + w] && reach[(size_t)w * n + v])
                redundant = 1;
        if (!redundant) ec[me++] = expr_copy((Expr*)e);
    }
    free(reach);
    graph_adj_free(a);

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vc, (size_t)n);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ec, me);
    free(vc); free(ec);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
