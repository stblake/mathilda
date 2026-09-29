/* components.c - connected-component builtins.
 *
 *   ConnectedComponents[g]          Mathematica's semantics, which depend on g:
 *                                     - g has a directed edge: the STRONGLY
 *                                       connected components, in an order with no
 *                                       edge from c_i to any later c_j (sinks
 *                                       first) -- Tarjan's completion order;
 *                                     - g is undirected: the components, largest
 *                                       first (ties by first appearance).
 *   ConnectedComponents[g, {v...}]  only the components containing some v.
 *   WeaklyConnectedComponents[g]    components of the underlying undirected
 *                                   graph, largest first; [g, {v...}] as above.
 *   StronglyConnectedComponents[g]  strong components (Tarjan) in the same
 *                                   sinks-first order as the directed
 *                                   ConnectedComponents (a Mathilda extension;
 *                                   Mathematica spells it ConnectedComponents).
 *
 * Every form lists the vertices of a component in VertexList order.
 *
 * Ordering, checked against Mathematica 15: Tarjan run from the vertices in
 * VertexList order, following out-edges in EdgeList order, emits the strong
 * components exactly in the order Mathematica returns them -- e.g.
 * ConnectedComponents[{3->1, 1->5, 2->4, 2->6, 3->5, 4->6}] is
 * {{5}, {1}, {3}, {6}, {4}, {2}} in both. Mathematica's order WITHIN an
 * undirected component, and between equal-sized undirected components, follows
 * no documented rule (it is not stable across such ties), so those use
 * VertexList order and first appearance respectively.
 *
 * Memory (SPEC section 4): results are freshly allocated; res is never
 * touched, so the evaluator frees it on success and keeps it on NULL.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Build the result List from a labeling comp[0..n-1] with `k` components,
 * listed in the order `order[0..k-1]` (component ids), keeping only those with
 * keep[cid] != 0 when keep is non-NULL. Vertices within a component keep
 * VertexList order. One counting pass, so O(n + k). */
static Expr* components_to_list(const GraphAdj* a, const int* comp, int k,
                                const int* order, const char* keep) {
    int n = a->n;
    size_t kk = (size_t)(k > 0 ? k : 1);
    int* cnt   = calloc(kk, sizeof(int));
    int* fill  = calloc(kk, sizeof(int));
    Expr*** members = calloc(kk, sizeof(Expr**));
    Expr** groups = calloc(kk, sizeof(Expr*));
    Expr* out = NULL;
    if (!cnt || !fill || !members || !groups) goto done;

    for (int i = 0; i < n; i++) cnt[comp[i]]++;
    for (int c = 0; c < k; c++) {
        if (keep && !keep[c]) continue;
        members[c] = malloc((size_t)(cnt[c] > 0 ? cnt[c] : 1) * sizeof(Expr*));
        if (!members[c]) goto done;
    }
    for (int i = 0; i < n; i++) {
        int c = comp[i];
        if (keep && !keep[c]) continue;
        members[c][fill[c]++] = expr_copy(a->verts->data.function.args[i]);
    }
    size_t ng = 0;
    for (int gi = 0; gi < k; gi++) {
        int c = order[gi];
        if (keep && !keep[c]) continue;
        groups[ng++] = expr_new_function(expr_new_symbol(SYM_List), members[c], (size_t)cnt[c]);
        free(members[c]);
        members[c] = NULL;             /* ownership moved into groups[] */
    }
    out = expr_new_function(expr_new_symbol(SYM_List), groups, ng);

done:
    if (members)
        for (int c = 0; c < k; c++)
            if (members[c]) {           /* only on the allocation-failure path */
                for (int j = 0; j < fill[c]; j++) expr_free(members[c][j]);
                free(members[c]);
            }
    free(cnt); free(fill); free(members); free(groups);
    return out;
}

/* Weak (underlying-undirected) labeling via DFS over out+in. Component ids are
 * assigned in order of each component's first vertex. Returns k, or -1 on
 * allocation failure. */
static int weak_label(const GraphAdj* a, int* comp) {
    int n = a->n;
    for (int i = 0; i < n; i++) comp[i] = -1;
    int* stack = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    if (!stack) return -1;
    int k = 0;
    for (int s = 0; s < n; s++) {
        if (comp[s] >= 0) continue;
        int top = 0; stack[top++] = s; comp[s] = k;
        while (top > 0) {
            int u = stack[--top];
            for (int j = 0; j < a->outdeg[u]; j++) { int w = a->out[u][j]; if (comp[w] < 0) { comp[w] = k; stack[top++] = w; } }
            for (int j = 0; j < a->indeg[u];  j++) { int w = a->in[u][j];  if (comp[w] < 0) { comp[w] = k; stack[top++] = w; } }
        }
        k++;
    }
    free(stack);
    return k;
}

/* ---- Tarjan strongly-connected components --------------------------------- */
typedef struct {
    const GraphAdj* a;
    int* index; int* low; char* onstack; int* stack; int sp;
    int* comp; int counter; int k;
} Tarjan;

/* Iterative Tarjan to avoid deep recursion on large graphs. Components are
 * numbered in completion order, which is a reverse topological order of the
 * condensation: no edge runs from component i to a component j > i. Returns 0
 * on allocation failure. */
static int tarjan_run(Tarjan* t) {
    int n = t->a->n;
    int* it = calloc((size_t)(n > 0 ? n : 1), sizeof(int)); /* per-node child cursor */
    int* callstk = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    if (!it || !callstk) { free(it); free(callstk); return 0; }
    for (int s = 0; s < n; s++) {
        if (t->index[s] >= 0) continue;
        int csp = 0; callstk[csp++] = s;
        t->index[s] = t->low[s] = t->counter++;
        t->stack[t->sp++] = s; t->onstack[s] = 1;
        while (csp > 0) {
            int u = callstk[csp - 1];
            if (it[u] < t->a->outdeg[u]) {
                int w = t->a->out[u][it[u]++];
                if (t->index[w] < 0) {
                    t->index[w] = t->low[w] = t->counter++;
                    t->stack[t->sp++] = w; t->onstack[w] = 1;
                    callstk[csp++] = w;
                } else if (t->onstack[w]) {
                    if (t->index[w] < t->low[u]) t->low[u] = t->index[w];
                }
            } else {
                /* Done with u: if root of an SCC, pop it. */
                if (t->low[u] == t->index[u]) {
                    int w;
                    do { w = t->stack[--t->sp]; t->onstack[w] = 0; t->comp[w] = t->k; } while (w != u);
                    t->k++;
                }
                csp--;
                if (csp > 0) { int p = callstk[csp - 1]; if (t->low[u] < t->low[p]) t->low[p] = t->low[u]; }
            }
        }
    }
    free(it); free(callstk);
    return 1;
}

int graph_strong_label(const GraphAdj* a, int* comp) {
    int n = a->n;
    size_t nn = (size_t)(n > 0 ? n : 1);
    Tarjan t;
    t.a = a;
    t.index = malloc(nn * sizeof(int));
    t.low   = malloc(nn * sizeof(int));
    t.onstack = calloc(nn, sizeof(char));
    t.stack = malloc(nn * sizeof(int));
    t.comp  = comp;
    t.sp = 0; t.counter = 0; t.k = 0;
    int ok = t.index && t.low && t.onstack && t.stack;
    if (ok) {
        for (int i = 0; i < n; i++) { t.index[i] = -1; t.low[i] = -1; comp[i] = -1; }
        ok = tarjan_run(&t);
    }
    free(t.index); free(t.low); free(t.onstack); free(t.stack);
    return ok ? t.k : -1;
}

/* ---- Shared driver --------------------------------------------------------- */

enum { CC_AUTO, CC_WEAK, CC_STRONG };

/* order[] for weak components: largest first, ties by first appearance (the
 * id order weak_label assigns). A stable counting sort on size, O(n + k). */
static int* order_by_size_desc(const int* comp, int n, int k) {
    size_t kk = (size_t)(k > 0 ? k : 1);
    int* size  = calloc(kk, sizeof(int));
    int* order = malloc(kk * sizeof(int));
    int* start = calloc((size_t)n + 2, sizeof(int));
    if (!size || !order || !start) { free(size); free(order); free(start); return NULL; }
    for (int i = 0; i < n; i++) size[comp[i]]++;
    /* bucket b = n - size, so bigger components come first */
    for (int c = 0; c < k; c++) start[n - size[c] + 1]++;
    for (int b = 1; b <= n + 1; b++) start[b] += start[b - 1];
    for (int c = 0; c < k; c++) order[start[n - size[c]]++] = c;
    free(size); free(start);
    return order;
}

static Expr* components_driver(Expr* res, int mode) {
    size_t argc = res->data.function.arg_count;
    if (argc != 1 && argc != 2) return NULL;
    const Expr* g = res->data.function.args[0];
    const Expr* sel = (argc == 2) ? res->data.function.args[1] : NULL;
    if (sel && !graph_is_list(sel)) return NULL;
    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    int n = a->n;

    int strong = (mode == CC_STRONG)
              || (mode == CC_AUTO && graph_directed_edge_count(g) > 0);
    int* comp = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    int* order = NULL;
    char* keep = NULL;
    Expr* out = NULL;
    int k = comp ? (strong ? graph_strong_label(a, comp) : weak_label(a, comp)) : -1;
    if (k < 0) goto done;

    if (strong) {
        /* Tarjan's completion order IS the answer's order. */
        order = malloc((size_t)(k > 0 ? k : 1) * sizeof(int));
        if (!order) goto done;
        for (int c = 0; c < k; c++) order[c] = c;
    } else {
        order = order_by_size_desc(comp, n, k);
        if (!order) goto done;
    }

    if (sel) {
        /* Keep the components containing at least one listed vertex; entries
         * that are not vertices of g select nothing. */
        keep = calloc((size_t)(k > 0 ? k : 1), 1);
        if (!keep) goto done;
        for (size_t i = 0; i < sel->data.function.arg_count; i++) {
            int p = graph_vertex_position(g, sel->data.function.args[i]);
            if (p >= 0) keep[comp[p]] = 1;
        }
    }
    out = components_to_list(a, comp, k, order, keep);

done:
    free(comp); free(order); free(keep);
    graph_adj_free(a);
    return out;
}

Expr* builtin_connected_components(Expr* res)        { return components_driver(res, CC_AUTO); }
Expr* builtin_weakly_connected_components(Expr* res) { return components_driver(res, CC_WEAK); }

Expr* builtin_strongly_connected_components(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    return components_driver(res, CC_STRONG);
}
