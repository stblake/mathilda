/* graphprops.c - structural property predicates.
 *
 *   UndirectedGraphQ[g]        every edge is undirected (edgeless graphs are
 *                              undirected, as in the Wolfram Language)
 *   EmptyGraphQ[g]             g has no edges (vertex count is irrelevant)
 *   CompleteGraphQ[g]          every ordered pair of distinct vertices (u, v) is
 *                              joined by an edge usable from u to v
 *   CompleteGraphQ[g, vlist]   the same test on the subgraph induced by vlist
 *   BipartiteGraphQ[g]         the underlying undirected graph is 2-colourable
 *
 * Like every *Q predicate, each gives False -- never unevaluated -- for an
 * argument that is not a valid graph.
 *
 * "Usable from u to v" follows GraphAdj.out[]: a DirectedEdge[u, v] or an
 * UndirectedEdge between u and v. A complete directed graph therefore needs
 * both u -> v and v -> u, matching Wolfram's "an undirected edge or two directed
 * edges going in opposite directions". Graphs with 0 or 1 vertices are complete
 * (there is no pair to fail).
 *
 * Memory (SPEC section 4): returns a fresh symbol; the evaluator frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

static Expr* truth(int b) { return expr_new_symbol(b ? SYM_True : SYM_False); }

Expr* builtin_undirected_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    return truth(graph_directed_edge_count(res->data.function.args[0]) == 0);
}

Expr* builtin_empty_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    return truth(graph_is_valid(g) && g->data.function.args[1]->data.function.arg_count == 0);
}

/* Completeness of the subgraph induced by the vertices with sel[i] != 0 (all of
 * them when sel is NULL). Each selected u must reach every other selected vertex
 * by an out-edge. out[u] may list a neighbour twice -- an UndirectedEdge u<->v
 * alongside a DirectedEdge u->v is legal -- so distinct successors are counted
 * with a per-u stamp rather than by degree. O(V + E). Returns 1/0, or -1 on
 * allocation failure (the caller then leaves the call unevaluated rather than
 * answer False). */
static int induced_complete(const GraphAdj* a, const char* sel) {
    int n = a->n, k = 0;
    for (int i = 0; i < n; i++) if (!sel || sel[i]) k++;
    if (k <= 1) return 1;
    int* stamp = malloc((size_t)n * sizeof(int));
    if (!stamp) return -1;
    for (int i = 0; i < n; i++) stamp[i] = -1;
    int ok = 1;
    for (int u = 0; u < n && ok; u++) {
        if (sel && !sel[u]) continue;
        int reach = 0;
        for (int j = 0; j < a->outdeg[u]; j++) {
            int w = a->out[u][j];
            if ((sel && !sel[w]) || stamp[w] == u) continue;
            stamp[w] = u;
            reach++;
        }
        if (reach != k - 1) ok = 0;
    }
    free(stamp);
    return ok;
}

Expr* builtin_complete_graph_q(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    long nd = graph_directed_edge_count(g);    /* O(1) via the memo; validates */
    if (nd < 0) return truth(0);

    if (argc == 1) {
        /* Simple graph with one edge kind: completeness is an edge count --
         * n(n-1)/2 undirected edges or n(n-1) directed ones, since parallel
         * edges are impossible. A mixed graph may pair u<->v with u->v, so it
         * takes the adjacency scan. */
        unsigned long long n = g->data.function.args[0]->data.function.arg_count;
        unsigned long long ne = g->data.function.args[1]->data.function.arg_count;
        unsigned long long und = (unsigned long long)nd;
        if (und == 0 || und == ne) {
            unsigned long long need = (n < 2) ? 0 : n * (n - 1) / (und == 0 ? 2 : 1);
            return truth(ne == need);
        }
        GraphAdj* a = graph_build_adj(g);
        if (!a) return NULL;
        int ok = induced_complete(a, NULL);
        graph_adj_free(a);
        return ok < 0 ? NULL : truth(ok);
    }

    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;

    /* CompleteGraphQ[g, vlist]: every listed item must be a vertex of g (else
     * False); repeats are harmless since sel[] is a set. */
    const Expr* vlist = res->data.function.args[1];
    if (!graph_is_list(vlist)) { graph_adj_free(a); return NULL; }
    char* sel = calloc((size_t)(a->n > 0 ? a->n : 1), sizeof(char));
    int ok = sel ? 1 : -1;
    for (size_t i = 0; ok > 0 && i < vlist->data.function.arg_count; i++) {
        int vi = graph_vertex_position(g, vlist->data.function.args[i]);
        if (vi < 0) ok = 0; else sel[vi] = 1;
    }
    if (ok > 0) ok = induced_complete(a, sel);
    free(sel);
    graph_adj_free(a);
    return ok < 0 ? NULL : truth(ok);
}

/* BFS 2-colouring over the underlying undirected graph (out[] and in[]), so
 * edge direction is ignored: bipartiteness is a property of the vertex split,
 * not of orientation. Edgeless graphs are trivially bipartite. */
Expr* builtin_bipartite_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return truth(0);
    int cached = graph_prop_get(g, GRAPH_PROP_BIPARTITE);   /* O(1) repeat */
    if (cached >= 0) return truth(cached);
    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    int n = a->n;
    signed char* side = malloc((size_t)(n > 0 ? n : 1));
    int* q = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!side || !q) { free(side); free(q); graph_adj_free(a); return NULL; }
    int ok = 1;
    for (int i = 0; i < n; i++) side[i] = -1;
    for (int s = 0; ok && s < n; s++) {
        if (side[s] >= 0) continue;
        int head = 0, tail = 0;
        side[s] = 0; q[tail++] = s;
        while (ok && head < tail) {
            int u = q[head++];
            for (int pass = 0; ok && pass < 2; pass++) {
                int deg = pass ? a->indeg[u] : a->outdeg[u];
                const int* nb = pass ? a->in[u] : a->out[u];
                for (int j = 0; j < deg; j++) {
                    int w = nb[j];
                    if (side[w] < 0) { side[w] = (signed char)(1 - side[u]); q[tail++] = w; }
                    else if (side[w] == side[u]) { ok = 0; break; }
                }
            }
        }
    }
    free(side); free(q);
    graph_adj_free(a);
    graph_prop_set(g, GRAPH_PROP_BIPARTITE, ok);
    return truth(ok);
}
