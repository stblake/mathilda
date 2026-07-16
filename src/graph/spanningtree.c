/* spanningtree.c - FindSpanningTree[g].
 *
 * A BFS spanning forest of the underlying undirected graph: for each component,
 * the tree edges chosen by BFS are collected in their original form (preserving
 * DirectedEdge/UndirectedEdge and orientation). Returns Graph[verts, treeEdges];
 * for a connected graph the tree has VertexCount - 1 edges.
 *
 * Memory (SPEC section 4): returns a freshly-allocated Graph; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Copy the original edge of g connecting vertex indices iu and iv (either
 * orientation / kind), or NULL if none. */
static Expr* original_edge_copy(const Expr* g, int iu, int iv) {
    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    for (size_t k = 0; k < edges->data.function.arg_count; k++) {
        const Expr* e = edges->data.function.args[k];
        int a = graph_vertex_index(verts, e->data.function.args[0]);
        int b = graph_vertex_index(verts, e->data.function.args[1]);
        if ((a == iu && b == iv) || (a == iv && b == iu)) return expr_copy((Expr*)e);
    }
    return NULL;
}

Expr* builtin_find_spanning_tree(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    int n = a->n;

    char* seen = calloc((size_t)(n > 0 ? n : 1), sizeof(char));
    int* q = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    Expr** tree = (n > 0) ? calloc((size_t)n, sizeof(Expr*)) : NULL;  /* <= n-1 edges */
    size_t te = 0;

    for (int s = 0; s < n; s++) {
        if (seen[s]) continue;
        int head = 0, tail = 0;
        seen[s] = 1; q[tail++] = s;
        while (head < tail) {
            int u = q[head++];
            /* Underlying undirected neighbors: out then in. */
            for (int pass = 0; pass < 2; pass++) {
                int deg = pass ? a->indeg[u] : a->outdeg[u];
                int* nb = pass ? a->in[u]   : a->out[u];
                for (int j = 0; j < deg; j++) {
                    int w = nb[j];
                    if (seen[w]) continue;
                    seen[w] = 1; q[tail++] = w;
                    Expr* e = original_edge_copy(g, u, w);
                    if (e) tree[te++] = e;
                }
            }
        }
    }

    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), tree, te);
    free(tree);
    /* Fresh copy of the vertex list. */
    size_t nv = (size_t)n;
    Expr** vcopy = (nv > 0) ? calloc(nv, sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < nv; i++) vcopy[i] = expr_copy(a->verts->data.function.args[i]);
    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), vcopy, nv);
    free(vcopy);

    free(seen); free(q); graph_adj_free(a);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
