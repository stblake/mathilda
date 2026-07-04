/* components.c - connected-component builtins.
 *
 *   ConnectedComponents[g]          weak components (underlying undirected)
 *   WeaklyConnectedComponents[g]    same as ConnectedComponents
 *   StronglyConnectedComponents[g]  strong components (Tarjan) over directed
 *                                   adjacency; for undirected graphs this
 *                                   coincides with the weak components.
 *
 * Each returns a List of Lists of vertices, components in first-appearance
 * order, vertices within a component in canonical index order.
 *
 * Memory (SPEC section 4): returns freshly-allocated lists; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* Build the result List from a component-id labeling comp[0..n-1] with `k`
 * components. Components appear in order of their first vertex; vertices within
 * each keep index order. */
static Expr* components_to_list(const GraphAdj* a, const int* comp, int k) {
    int n = a->n;
    /* First-appearance order of component ids. */
    int* order = calloc((size_t)(k > 0 ? k : 1), sizeof(int));
    int* pos   = malloc((size_t)(k > 0 ? k : 1) * sizeof(int));
    for (int i = 0; i < k; i++) pos[i] = -1;
    int seen = 0;
    for (int i = 0; i < n; i++)
        if (pos[comp[i]] < 0) { pos[comp[i]] = seen; order[seen++] = comp[i]; }

    Expr** groups = calloc((size_t)(k > 0 ? k : 1), sizeof(Expr*));
    for (int gi = 0; gi < k; gi++) {
        int cid = order[gi];
        int cnt = 0;
        for (int i = 0; i < n; i++) if (comp[i] == cid) cnt++;
        Expr** members = (cnt > 0) ? calloc((size_t)cnt, sizeof(Expr*)) : NULL;
        int m = 0;
        for (int i = 0; i < n; i++)
            if (comp[i] == cid) members[m++] = expr_copy(a->verts->data.function.args[i]);
        groups[gi] = expr_new_function(expr_new_symbol(SYM_List), members, (size_t)cnt);
        free(members);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), groups, (size_t)k);
    free(groups); free(order); free(pos);
    return out;
}

/* Weak/underlying-undirected labeling via DFS over out+in. */
static Expr* weak_components(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GraphAdj* a = graph_build_adj(res->data.function.args[0]);
    if (!a) return NULL;
    int n = a->n;
    int* comp = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    for (int i = 0; i < n; i++) comp[i] = -1;
    int* stack = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
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
    Expr* out = components_to_list(a, comp, k);
    free(comp); free(stack); graph_adj_free(a);
    return out;
}

Expr* builtin_connected_components(Expr* res)          { return weak_components(res); }
Expr* builtin_weakly_connected_components(Expr* res)   { return weak_components(res); }

/* ---- Tarjan strongly-connected components --------------------------------- */
typedef struct {
    const GraphAdj* a;
    int* index; int* low; char* onstack; int* stack; int sp;
    int* comp; int counter; int k;
} Tarjan;

/* Iterative Tarjan to avoid deep recursion on large graphs. */
static void tarjan_run(Tarjan* t) {
    int n = t->a->n;
    int* it = calloc((size_t)(n > 0 ? n : 1), sizeof(int)); /* per-node child cursor */
    int* callstk = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
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
}

Expr* builtin_strongly_connected_components(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    GraphAdj* a = graph_build_adj(res->data.function.args[0]);
    if (!a) return NULL;
    int n = a->n;
    Tarjan t;
    t.a = a;
    t.index = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    t.low   = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    t.onstack = calloc((size_t)(n > 0 ? n : 1), sizeof(char));
    t.stack = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    t.comp  = malloc((size_t)(n > 0 ? n : 1) * sizeof(int));
    t.sp = 0; t.counter = 0; t.k = 0;
    for (int i = 0; i < n; i++) { t.index[i] = -1; t.low[i] = -1; t.comp[i] = -1; }
    tarjan_run(&t);
    Expr* out = components_to_list(a, t.comp, t.k);
    free(t.index); free(t.low); free(t.onstack); free(t.stack); free(t.comp);
    graph_adj_free(a);
    return out;
}
