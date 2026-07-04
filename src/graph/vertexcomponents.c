/* vertexcomponents.c - VertexOutComponent[g, v] and VertexInComponent[g, v].
 *
 *   VertexOutComponent[g, v]  the vertices reachable FROM v following edges
 *                             forward (over out-edges), including v itself.
 *   VertexInComponent[g, v]   the vertices from which v is reachable (over
 *                             in-edges), including v itself.
 *
 * Each is a single BFS from v over the appropriate direction of the adjacency,
 * O(V+E); results are vertices in the graph's canonical order. For an undirected
 * graph out- and in-neighbourhoods coincide, so both give v's connected
 * component. A vertex not in g leaves the call unevaluated.
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph
 * or an unknown vertex.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Shared: BFS from vertex `res`'s second argument over out[] (use_in=0) or
 * in[] (use_in=1); returns the reachable vertices as a List, or NULL. */
static Expr* component(Expr* res, int use_in) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    const Expr* verts = g->data.function.args[0];
    int src = graph_vertex_index(verts, res->data.function.args[1]);
    if (src < 0) return NULL;

    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    int n = a->n;
    int**  nbr = use_in ? a->in     : a->out;
    int*   deg = use_in ? a->indeg  : a->outdeg;

    char* seen = calloc((size_t)n, 1);
    int*  q    = malloc((size_t)n * sizeof(int));
    int head = 0, tail = 0;
    seen[src] = 1; q[tail++] = src;
    while (head < tail) {
        int v = q[head++];
        for (int j = 0; j < deg[v]; j++) {
            int u = nbr[v][j];
            if (!seen[u]) { seen[u] = 1; q[tail++] = u; }
        }
    }
    free(q);
    graph_adj_free(a);

    Expr** items = malloc((size_t)n * sizeof(Expr*));
    size_t cnt = 0;
    for (int i = 0; i < n; i++) if (seen[i]) items[cnt++] = expr_copy(verts->data.function.args[i]);
    free(seen);
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, cnt);
    free(items);
    return out;
}

Expr* builtin_vertex_out_component(Expr* res) { return component(res, 0); }
Expr* builtin_vertex_in_component(Expr* res)  { return component(res, 1); }
