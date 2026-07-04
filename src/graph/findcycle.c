/* findcycle.c - FindCycle[g]: a cycle in g as a list containing one cycle, that
 * cycle being a list of its edges (Wolfram's shape, e.g. {{1<->2, 2<->3,
 * 3<->1}}), or {} when g is acyclic.
 *
 * Depth-first search with back-edge detection, O(V+E):
 *   directed g   -> a gray (on-stack) target closes a directed cycle;
 *   undirected g -> a visited neighbour that is not the DFS parent closes one
 *                   (undirected DFS has only tree and back edges, so any such
 *                   neighbour is a genuine ancestor).
 * The cycle is rebuilt by walking DFS-tree parent pointers from the back edge's
 * tail up to its head. The first cycle found (deterministic: sources and
 * neighbours in canonical order) is returned -- not necessarily the shortest.
 *
 * Edge kinds mirror the source graph: DirectedEdge for a fully directed graph,
 * UndirectedEdge otherwise. A cycle needs >= 3 vertices in the undirected case,
 * automatic under the simple-graph invariant (no self-loops / parallel edges).
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

typedef struct {
    const GraphAdj* a;
    int*  state;    /* directed: 0 white / 1 gray / 2 black                    */
    char* visited;  /* undirected                                              */
    int*  parent;   /* DFS-tree parent (vertex index), -1 at a root            */
    int*  cyc;      /* reconstructed cycle vertices, length cyclen             */
    int*  tmp;      /* scratch for reconstruction                              */
    int   cyclen;   /* 0 while no cycle found                                  */
} CycEnv;

/* Record the cycle head..tail (ancestor u down to back-edge tail v). */
static void cyc_record(CycEnv* e, int v, int u) {
    int len = 0, w = v;
    for (;;) { e->tmp[len++] = w; if (w == u) break; w = e->parent[w]; }
    for (int i = 0; i < len; i++) e->cyc[i] = e->tmp[len - 1 - i]; /* u..v */
    e->cyclen = len;
}

static int dfs_dir(CycEnv* e, int v) {
    e->state[v] = 1;
    for (int j = 0; j < e->a->outdeg[v]; j++) {
        int u = e->a->out[v][j];
        if (e->state[u] == 0) { e->parent[u] = v; if (dfs_dir(e, u)) return 1; }
        else if (e->state[u] == 1) { cyc_record(e, v, u); return 1; }
    }
    e->state[v] = 2;
    return 0;
}

static int dfs_undir(CycEnv* e, int v, int par) {
    e->visited[v] = 1; e->parent[v] = par;
    for (int j = 0; j < e->a->outdeg[v]; j++) {
        int u = e->a->out[v][j];
        if (!e->visited[u]) { if (dfs_undir(e, u, v)) return 1; }
        else if (u != par) { cyc_record(e, v, u); return 1; }
    }
    return 0;
}

Expr* builtin_find_cycle(Expr* res) {
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

    CycEnv e;
    e.a = a;
    e.state   = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    e.visited = calloc((size_t)(n > 0 ? n : 1), 1);
    e.parent  = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    e.cyc     = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    e.tmp     = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    e.cyclen  = 0;
    for (int i = 0; i < n; i++) e.parent[i] = -1;

    for (int s = 0; s < n && e.cyclen == 0; s++) {
        if (directed) { if (e.state[s] == 0) dfs_dir(&e, s); }
        else          { if (!e.visited[s])   dfs_undir(&e, s, -1); }
    }

    Expr* out;
    if (e.cyclen > 0) {
        const char* ekind = directed ? SYM_DirectedEdge : SYM_UndirectedEdge;
        Expr** cedges = calloc((size_t)e.cyclen, sizeof(Expr*));
        for (int i = 0; i < e.cyclen; i++) {
            int from = e.cyc[i], to = e.cyc[(i + 1) % e.cyclen];
            Expr* args[2] = { expr_copy(verts->data.function.args[from]),
                              expr_copy(verts->data.function.args[to]) };
            cedges[i] = expr_new_function(expr_new_symbol(ekind), args, 2);
        }
        Expr* cycle = expr_new_function(expr_new_symbol(SYM_List), cedges, (size_t)e.cyclen);
        free(cedges);
        Expr* one[1] = { cycle };
        out = expr_new_function(expr_new_symbol(SYM_List), one, 1);
    } else {
        out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    free(e.state); free(e.visited); free(e.parent); free(e.cyc); free(e.tmp);
    graph_adj_free(a);
    return out;
}
