/* transitiveclosure.c - TransitiveClosure[g].
 *
 *   directed g   -> directed closure: an edge u->v for every v (!= u) reachable
 *                   from u by a directed path.
 *   undirected g -> each connected component becomes a clique (u<->v for every
 *                   pair in the same component).
 *
 * Directed: a BFS per source over out[] marks its reachable set (O(V*(V+E))).
 * Undirected: one component labeling, then all within-component pairs.
 *
 * Memory (SPEC section 4): returns the canonical Graph directly; frees res.
 * NULL on a non-graph argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_transitive_closure(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;
    int directed = (ne > 0);
    for (size_t i = 0; i < ne; i++)
        if (graph_edge_kind(edges->data.function.args[i]) != SYM_DirectedEdge) { directed = 0; break; }

    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;

    /* Worst case n(n-1) directed edges (n(n-1)/2 undirected). */
    size_t cap = (size_t)n * (size_t)(n > 0 ? n - 1 : 0);
    Expr** cedges = (cap > 0) ? calloc(cap, sizeof(Expr*)) : NULL;
    size_t m = 0;
    int* dist = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;
    int* q    = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;

    if (directed) {
        for (int s = 0; s < n; s++) {
            for (int i = 0; i < n; i++) dist[i] = -1;
            int head = 0, tail = 0; dist[s] = 0; q[tail++] = s;
            while (head < tail) {
                int u = q[head++];
                for (int e = 0; e < a->outdeg[u]; e++) {
                    int w = a->out[u][e];
                    if (dist[w] < 0) { dist[w] = dist[u] + 1; q[tail++] = w; }
                }
            }
            for (int v = 0; v < n; v++)
                if (v != s && dist[v] >= 0) {
                    Expr* ea[2] = { expr_copy((Expr*)verts->data.function.args[s]),
                                    expr_copy((Expr*)verts->data.function.args[v]) };
                    cedges[m++] = expr_new_function(expr_new_symbol(SYM_DirectedEdge), ea, 2);
                }
        }
    } else {
        /* Undirected: label components (BFS over out[] + in[]), then clique each. */
        int* comp = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;
        for (int i = 0; i < n; i++) comp[i] = -1;
        int nc = 0;
        for (int s = 0; s < n; s++) {
            if (comp[s] >= 0) continue;
            int head = 0, tail = 0; comp[s] = nc; q[tail++] = s;
            while (head < tail) {
                int u = q[head++];
                for (int pass = 0; pass < 2; pass++) {
                    int deg = pass ? a->indeg[u] : a->outdeg[u];
                    int* nb = pass ? a->in[u] : a->out[u];
                    for (int e = 0; e < deg; e++)
                        if (comp[nb[e]] < 0) { comp[nb[e]] = nc; q[tail++] = nb[e]; }
                }
            }
            nc++;
        }
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (comp[i] == comp[j]) {
                    Expr* ea[2] = { expr_copy((Expr*)verts->data.function.args[i]),
                                    expr_copy((Expr*)verts->data.function.args[j]) };
                    cedges[m++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), ea, 2);
                }
        free(comp);
    }
    free(dist); free(q);

    Expr** vcopy = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    for (int i = 0; i < n; i++) vcopy[i] = expr_copy((Expr*)verts->data.function.args[i]);
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vcopy, (size_t)n);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), cedges, m);
    free(vcopy); free(cedges);
    graph_adj_free(a);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
