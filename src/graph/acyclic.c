/* acyclic.c - cycle structure and ordering.
 *
 *   AcyclicGraphQ[g]           True iff g has no cycle, respecting direction
 *   TreeGraphQ[g]              True iff g is a tree: connected, >= 1 vertex,
 *                              and exactly VertexCount - 1 edges
 *   TopologicalSort[g]         vertices of a directed acyclic graph with u
 *                              before v for every edge u -> v
 *   TopologicalSort[{rules}]   the same, for Graph[{rules}]
 *
 * AcyclicGraphQ -- a cycle is a closed walk that uses no edge twice, following
 * directed edges forwards only and undirected edges either way. So an
 * undirected graph is acyclic iff it is a forest, a directed graph iff it is a
 * DAG, and u -> v with v -> u is a (2-)cycle. Mixed graphs are decided exactly
 * by contraction:
 *   1. Union-find over the undirected edges. Joining two vertices already in
 *      one set closes an undirected cycle.
 *   2. Otherwise every undirected component is a tree, so any vertex in it can
 *      reach any other without reusing an edge. Collapse each to one node. A
 *      directed edge inside a component closes a cycle with the tree path back
 *      to its tail; the rest form a directed multigraph on the components.
 *   3. g is acyclic iff that contracted digraph is (Kahn's algorithm). A simple
 *      cycle there visits each component once, so it lifts to a cycle of g
 *      using one tree path per component; conversely any cycle of g that uses
 *      a directed edge projects to a closed walk in the contraction.
 * Linear time, O(V + E alpha(V)).
 *
 * All three cache their answer on the graph's validated-graph memo entry
 * (graph_prop_*, graph_cached_*), so repeating a query on the same graph node
 * is O(1), as with Mathematica's atomic Graph object.
 *
 * TreeGraphQ follows Wolfram's "a simple connected graph with no cycles" on the
 * underlying undirected graph, where "connected with n - 1 edges" is the
 * equivalent test (connectivity by union-find over the edges). Direction is ignored, so an out-tree like 1->2, 1->3 is a
 * tree; the anti-parallel pair 1->2, 2->1 is two edges on two vertices, not a
 * tree. The null graph (no vertices) is not a tree, matching ConnectedGraphQ.
 *
 * TopologicalSort is Kahn's algorithm with ties broken by VertexList position
 * (the smallest ready index goes first, via a binary min-heap), so the order is
 * deterministic and an already-sorted vertex list is returned unchanged. It is
 * left unevaluated for a graph with a cycle, or with any undirected edge (an
 * undirected edge imposes no order, so the graph is not a DAG); an edgeless
 * graph sorts to its VertexList.
 *
 * AcyclicGraphQ/TreeGraphQ give False for a non-graph; TopologicalSort is left
 * unevaluated.
 *
 * Memory (SPEC section 4): returns freshly-allocated results; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include <stdlib.h>

static Expr* truth(int b) { return expr_new_symbol(b ? SYM_True : SYM_False); }

/* ---- Union-find (path halving + union by size) ---------------------------- */
static int uf_find(int* parent, int x) {
    while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
    return x;
}

/* Joins the sets of a and b; returns 0 if they were already one set. */
static int uf_union(int* parent, int* size, int a, int b) {
    a = uf_find(parent, a); b = uf_find(parent, b);
    if (a == b) return 0;
    if (size[a] < size[b]) { int t = a; a = b; b = t; }
    parent[b] = a; size[a] += size[b];
    return 1;
}

/* 1 if the valid graph g is acyclic, 0 if it has a cycle, -1 if g is not a
 * valid graph or an allocation failed. Endpoints come pre-resolved from the
 * validated-graph memo, so this is integer-only work. */
static int graph_acyclic(const Expr* g) {
    const int *eu, *ev;
    const unsigned char* edir;
    if (!graph_edge_indices(g, &eu, &ev, &edir)) return -1;
    int n = (int)g->data.function.args[0]->data.function.arg_count;
    size_t ne = g->data.function.args[1]->data.function.arg_count;
    size_t nn = (size_t)(n > 0 ? n : 1);

    int* parent = malloc(nn * sizeof(int));
    int* size   = malloc(nn * sizeof(int));
    int* ea     = malloc((ne > 0 ? ne : 1) * sizeof(int));   /* directed tails */
    int* eb     = malloc((ne > 0 ? ne : 1) * sizeof(int));   /* directed heads */
    int* indeg  = calloc(nn, sizeof(int));
    int* outdeg = calloc(nn, sizeof(int));
    int result = -1;
    if (!parent || !size || !ea || !eb || !indeg || !outdeg) goto done;

    for (int i = 0; i < n; i++) { parent[i] = i; size[i] = 1; }

    /* Step 1: undirected edges -- any repeat join is a cycle. Directed edges
     * are staged for step 2, which needs the final components. */
    size_t nd = 0;
    result = 1;
    for (size_t k = 0; k < ne && result; k++) {
        if (!edir[k]) {
            if (!uf_union(parent, size, eu[k], ev[k])) result = 0;
        } else {
            ea[nd] = eu[k]; eb[nd] = ev[k]; nd++;
        }
    }
    if (!result) goto done;

    /* Step 2: map directed edges onto components; an intra-component edge
     * closes a cycle through the component's tree. */
    for (size_t k = 0; k < nd; k++) {
        ea[k] = uf_find(parent, ea[k]);
        eb[k] = uf_find(parent, eb[k]);
        if (ea[k] == eb[k]) { result = 0; goto done; }
        outdeg[ea[k]]++; indeg[eb[k]]++;
    }

    /* Step 3: Kahn over the contracted digraph (nodes = component roots), in
     * CSR form. Acyclic iff every root is eventually removed. */
    {
        int* start = calloc(nn + 1, sizeof(int));
        int* succ  = malloc((nd > 0 ? nd : 1) * sizeof(int));
        int* queue = malloc(nn * sizeof(int));
        if (!start || !succ || !queue) {
            free(start); free(succ); free(queue); result = -1; goto done;
        }
        for (int i = 0; i < n; i++) start[i + 1] = start[i] + outdeg[i];
        for (size_t k = 0; k < nd; k++) succ[start[ea[k]] + --outdeg[ea[k]]] = eb[k];
        int roots = 0, head = 0, tail = 0;
        for (int i = 0; i < n; i++) {
            if (uf_find(parent, i) != i) continue;
            roots++;
            if (indeg[i] == 0) queue[tail++] = i;
        }
        while (head < tail) {
            int u = queue[head++];
            for (int j = start[u]; j < start[u + 1]; j++)
                if (--indeg[succ[j]] == 0) queue[tail++] = succ[j];
        }
        result = (tail == roots);
        free(start); free(succ); free(queue);
    }

done:
    free(parent); free(size); free(ea); free(eb); free(indeg); free(outdeg);
    return result;
}

Expr* builtin_acyclic_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return truth(0);
    int r = graph_prop_get(g, GRAPH_PROP_ACYCLIC);
    if (r < 0) {
        r = graph_acyclic(g);
        if (r < 0) return NULL;                      /* allocation failure */
        graph_prop_set(g, GRAPH_PROP_ACYCLIC, r);
    }
    return truth(r);
}

/* A tree iff >= 1 vertex, n - 1 edges and connected. Given the edge count, the
 * connectivity test is a union-find over the edge indices -- exactly n - 1
 * successful joins -- which needs no adjacency at all. */
Expr* builtin_tree_graph_q(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    const int *eu, *ev;
    const unsigned char* edir;
    if (!graph_edge_indices(g, &eu, &ev, &edir)) return truth(0);
    int cached = graph_prop_get(g, GRAPH_PROP_TREE);
    if (cached >= 0) return truth(cached);

    int n = (int)g->data.function.args[0]->data.function.arg_count;
    size_t ne = g->data.function.args[1]->data.function.arg_count;
    int tree = 0;
    if (n >= 1 && ne == (size_t)(n - 1)) {
        int* parent = malloc((size_t)n * sizeof(int));
        int* size   = malloc((size_t)n * sizeof(int));
        if (!parent || !size) { free(parent); free(size); return NULL; }
        for (int i = 0; i < n; i++) { parent[i] = i; size[i] = 1; }
        size_t joins = 0;
        for (size_t k = 0; k < ne; k++) joins += (size_t)uf_union(parent, size, eu[k], ev[k]);
        tree = (joins == ne);          /* n - 1 joins leave one component */
        free(parent); free(size);
    }
    graph_prop_set(g, GRAPH_PROP_TREE, tree);
    return truth(tree);
}

/* ---- TopologicalSort ------------------------------------------------------ */

/* Binary min-heap of vertex indices, for smallest-position-first tie-breaking. */
static void heap_push(int* h, int* len, int x) {
    int i = (*len)++;
    while (i > 0) {
        int p = (i - 1) / 2;
        if (h[p] <= x) break;
        h[i] = h[p]; i = p;
    }
    h[i] = x;
}

static int heap_pop(int* h, int* len) {
    int top = h[0], x = h[--(*len)], i = 0;
    for (;;) {
        int c = 2 * i + 1;
        if (c >= *len) break;
        if (c + 1 < *len && h[c + 1] < h[c]) c++;
        if (x <= h[c]) break;
        h[i] = h[c]; i = c;
    }
    if (*len > 0) h[i] = x;
    return top;
}

/* Kahn over a validated graph. Returns the sorted vertex List, or NULL when g
 * has an undirected edge, has a cycle, or an allocation failed. */
static Expr* topo_sort(const Expr* g) {
    const Expr* edges = g->data.function.args[1];
    for (size_t k = 0; k < edges->data.function.arg_count; k++)
        if (graph_edge_kind(edges->data.function.args[k]) != SYM_DirectedEdge)
            return NULL;

    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    int n = a->n;
    size_t nn = (size_t)(n > 0 ? n : 1);
    int* indeg = malloc(nn * sizeof(int));
    int* heap  = malloc(nn * sizeof(int));
    Expr** out = calloc(nn, sizeof(Expr*));
    Expr* result = NULL;
    if (!indeg || !heap || !out) goto done;

    int len = 0, done_count = 0;
    for (int i = 0; i < n; i++) {
        indeg[i] = a->indeg[i];
        if (indeg[i] == 0) heap_push(heap, &len, i);
    }
    while (len > 0) {
        int u = heap_pop(heap, &len);
        out[done_count++] = expr_copy(a->verts->data.function.args[u]);
        for (int j = 0; j < a->outdeg[u]; j++)
            if (--indeg[a->out[u][j]] == 0) heap_push(heap, &len, a->out[u][j]);
    }
    if (done_count == n) {
        result = expr_new_function(expr_new_symbol(SYM_List), out, (size_t)n);
    } else {
        for (int i = 0; i < done_count; i++) expr_free(out[i]);   /* cyclic */
    }

done:
    free(indeg); free(heap); free(out);
    graph_adj_free(a);
    return result;
}

/* topo_sort with the result cached on g's memo entry (a list is an immutable
 * value, so handing out references to one is safe). A failure (cycle or
 * undirected edge) is recorded as "not acyclic-and-directed" by recomputing --
 * rare, and still linear. */
static Expr* topo_sort_cached(const Expr* g) {
    Expr* hit = graph_cached_get(g, GRAPH_CACHED_TOPOSORT);
    if (hit) return hit;
    Expr* out = topo_sort(g);
    if (out) graph_cached_set(g, GRAPH_CACHED_TOPOSORT, out);
    return out;
}

Expr* builtin_topological_sort(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* arg = res->data.function.args[0];
    if (graph_is_valid(arg)) return topo_sort_cached(arg);
    if (!graph_is_list(arg)) return NULL;

    /* TopologicalSort[{rules}]: build Graph[{rules}] through the constructor,
     * so the edge sugar and validation are exactly Graph's. */
    Expr* garg = expr_copy(res->data.function.args[0]);
    Expr* call = expr_new_function(expr_new_symbol(SYM_Graph), &garg, 1);
    Expr* g = evaluate(call);            /* evaluate() does not consume call */
    expr_free(call);
    Expr* out = graph_is_valid(g) ? topo_sort(g) : NULL;
    expr_free(g);
    return out;
}
