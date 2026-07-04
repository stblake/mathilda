/* graphpower.c - GraphPower[g, k]: the k-th power of a graph. Two vertices are
 * joined in the result iff they are connected by a path of length 1..k in g
 * (a vertex is never joined to itself -- no self-loops).
 *
 *   undirected g -> undirected power: {i,j} joined iff dist(i,j) <= k.
 *   directed g   -> directed power:   i->j    iff the directed dist(i,j) <= k.
 *
 * Reachability comes from a depth-limited BFS per source over GraphAdj.out[],
 * which already encodes direction (an undirected edge sits in both endpoints'
 * out[]), so one traversal handles both cases. Cost O(V*(V+E)); the result is
 * the canonical Graph value, ready to render or query.
 *
 * k must be a positive integer literal; otherwise the call is left unevaluated
 * (NULL) so GraphPower[g, k] stays symbolic for non-numeric k.
 *
 * Memory (SPEC section 4): returns a freshly-built Graph tree; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_graph_power(Expr* res) {
    if (res->data.function.arg_count != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    const Expr* ke = res->data.function.args[1];
    if (!graph_is_valid(g)) return NULL;
    if (ke->type != EXPR_INTEGER || ke->data.integer < 1) return NULL;
    long k = (long)ke->data.integer;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    int directed = (ne > 0);
    for (size_t i = 0; i < ne; i++)
        if (graph_edge_kind(edges->data.function.args[i]) != SYM_DirectedEdge) { directed = 0; break; }

    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;

    /* n(n-1) is the dense upper bound (undirected halves it). */
    size_t cap = (size_t)n * (n > 0 ? (size_t)(n - 1) : 0);
    if (!directed) cap /= 2;
    Expr** pedges = (cap > 0) ? calloc(cap, sizeof(Expr*)) : NULL;
    size_t m = 0;
    const char* ekind = directed ? SYM_DirectedEdge : SYM_UndirectedEdge;

    int* dist  = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;
    int* queue = (n > 0) ? malloc((size_t)n * sizeof(int)) : NULL;
    for (int s = 0; s < n; s++) {
        for (int i = 0; i < n; i++) dist[i] = -1;
        int qh = 0, qt = 0;
        dist[s] = 0; queue[qt++] = s;
        while (qh < qt) {
            int v = queue[qh++];
            if (dist[v] >= (int)k) continue;         /* depth-limited */
            for (int j = 0; j < a->outdeg[v]; j++) {
                int u = a->out[v][j];
                if (dist[u] < 0) { dist[u] = dist[v] + 1; queue[qt++] = u; }
            }
        }
        for (int t = 0; t < n; t++) {
            if (t == s || dist[t] < 0) continue;     /* unreachable / self */
            if (!directed && t < s) continue;        /* undirected: emit once */
            Expr* args[2] = { expr_copy(verts->data.function.args[s]),
                              expr_copy(verts->data.function.args[t]) };
            pedges[m++] = expr_new_function(expr_new_symbol(ekind), args, 2);
        }
    }
    free(dist); free(queue);
    graph_adj_free(a);

    Expr** vcopy = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;
    for (int i = 0; i < n; i++) vcopy[i] = expr_copy(verts->data.function.args[i]);
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vcopy, (size_t)n);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), pedges, m);
    free(vcopy); free(pedges);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
