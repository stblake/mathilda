/* gops_preds.c - structural predicates.
 *
 *   SimpleGraphQ[g]        True for every valid graph (Mathilda graphs have no
 *                          self-loops or parallel edges); False otherwise.
 *   LoopFreeGraphQ[g]      likewise True for every valid graph.
 *   MixedGraphQ[g]         True iff g has both directed and undirected edges.
 *   WeightedGraphQ[g]      True iff g carries EdgeWeight (Mathilda has no
 *   EdgeWeightedGraphQ[g]  vertex weights, so the two agree).
 *   PathGraphQ[g]          Mathematica's definition: at least one vertex,
 *                          connected, and
 *                            undirected: every degree <= 2;
 *                            directed:   every in- and out-degree <= 1;
 *                          so a cycle counts as a (closed) path, as it does in
 *                          Mathematica 15 (PathGraphQ[CycleGraph[3]] is True);
 *                            mixed:      never.
 *   EulerianGraphQ[g]      True iff g has an Eulerian cycle: every vertex of
 *                          even degree (undirected) or with in-degree =
 *                          out-degree (directed), and all edges in one
 *                          connected component. Edgeless graphs with >= 1
 *                          vertex are Eulerian; the null graph is not. Mixed
 *                          graphs are left unevaluated (as Mathematica).
 *
 * The *Q heads give False for a non-graph (EulerianGraphQ included). All are
 * O(1) or one O(V + E) integer pass over the memo's endpoints.
 *
 * Memory (SPEC section 4): results are fresh symbols; res is borrowed.
 */

#include "graph_ops.h"
#include <stdlib.h>
#include <string.h>

Expr* builtin_simple_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    return gops_truth(graph_is_valid(res->data.function.args[0]));
}

Expr* builtin_loop_free_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    return gops_truth(graph_is_valid(res->data.function.args[0]));
}

Expr* builtin_mixed_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    long nd = graph_directed_edge_count(g);
    if (nd < 0) return gops_truth(0);
    size_t ne = g->data.function.args[1]->data.function.arg_count;
    return gops_truth(nd > 0 && (size_t)nd < ne);
}

Expr* builtin_weighted_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    return gops_truth(graph_is_valid(g) && g->data.function.arg_count == 3);
}

Expr* builtin_edge_weighted_graph_q(Expr* res) {
    return builtin_weighted_graph_q(res);
}

/* Union-find with path halving. */
static int uf_find(int* p, int x) {
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}

/* Number of union-find joins over all edges (so the graph restricted to its
 * non-isolated vertices is connected iff joins == touched - 1). */
static size_t count_joins(const GopsView* v, int* parent) {
    size_t joins = 0;
    for (size_t i = 0; i < v->nv; i++) parent[i] = (int)i;
    for (size_t k = 0; k < v->ne; k++) {
        int a = uf_find(parent, v->eu[k]), b = uf_find(parent, v->ev[k]);
        if (a != b) { parent[a] = b; joins++; }
    }
    return joins;
}

Expr* builtin_path_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GopsView v;
    if (!gops_view(res->data.function.args[0], &v)) return gops_truth(0);
    size_t n = v.nv, ne = v.ne;
    if (n == 0) return gops_truth(0);
    if (v.ndir != 0 && v.ndir != ne) return gops_truth(0);     /* mixed */
    int directed = (v.ndir == ne && ne > 0);
    if (ne > n) return gops_truth(0);         /* degree bounds allow <= n edges */
    int* din = gops_calloc(n, sizeof(int));
    int* dout = gops_calloc(n, sizeof(int));
    int* parent = gops_malloc(n, sizeof(int));
    if (!din || !dout || !parent) { free(din); free(dout); free(parent); return NULL; }
    int ok = 1;
    for (size_t k = 0; k < ne && ok; k++) {
        int a = v.eu[k], b = v.ev[k];
        if (directed) {
            if (++dout[a] > 1 || ++din[b] > 1) ok = 0;
        } else {
            if (++dout[a] + din[a] > 2 || ++din[b] + dout[b] > 2) ok = 0;
        }
    }
    if (ok) ok = (count_joins(&v, parent) == n - 1);           /* connected */
    free(din); free(dout); free(parent);
    return gops_truth(ok);
}

/* 1 Eulerian, 0 not, -1 mixed (unevaluated), -2 allocation failure. */
int gops_eulerian(const GopsView* v) {
    size_t n = v->nv, ne = v->ne;
    if (n == 0) return 0;
    if (ne == 0) return 1;
    if (v->ndir != 0 && v->ndir != ne) return -1;
    int directed = (v->ndir == ne);
    int* bal = gops_calloc(n, sizeof(int));       /* degree parity / balance */
    unsigned char* touched = gops_calloc(n, 1);
    int* parent = gops_malloc(n, sizeof(int));
    if (!bal || !touched || !parent) { free(bal); free(touched); free(parent); return -2; }
    for (size_t k = 0; k < ne; k++) {
        int a = v->eu[k], b = v->ev[k];
        touched[a] = touched[b] = 1;
        if (directed) { bal[a]++; bal[b]--; }
        else { bal[a] ^= 1; bal[b] ^= 1; }
    }
    int ok = 1;
    size_t nt = 0;
    for (size_t i = 0; i < n && ok; i++) { if (bal[i]) ok = 0; nt += touched[i]; }
    /* Balanced + weakly connected implies strongly connected (directed). */
    if (ok) ok = (count_joins(v, parent) == nt - 1);
    free(bal); free(touched); free(parent);
    return ok;
}

Expr* builtin_eulerian_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GopsView v;
    if (!gops_view(res->data.function.args[0], &v)) return gops_truth(0);
    int r = gops_eulerian(&v);
    if (r < 0) return NULL;
    return gops_truth(r);
}
