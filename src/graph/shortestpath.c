/* shortestpath.c - FindShortestPath[g,s,t] and GraphDistance[g,s,t].
 *
 * Two algorithms, dispatched on graph_weights_usable(g):
 *   - Unweighted (default): breadth-first search over the successor adjacency
 *     (GraphAdj.out[]): for a directed graph this follows edge direction; for
 *     an undirected graph out[] is symmetric, so it is an ordinary shortest
 *     path.
 *   - Weighted (g carries a non-negative-numeric EdgeWeight): Dijkstra over a
 *     local, call-scoped weighted adjacency (WAdj below) -- NOT over
 *     GraphAdj, which has no weight storage and is shared by 5 other
 *     builtins (ConnectedComponents, WeaklyConnectedComponents,
 *     FindSpanningTree, ConnectedGraphQ, VertexConnectivity); widening it
 *     would risk the exact class of shared-choke-point defect a
 *     plan-reviewer pass caught during the EdgeWeight ticket. Falls back to
 *     BFS for a symbolic or negative weight rather than erroring.
 *
 * Wolfram's naming split is kept: FindShortestPath returns the vertex path,
 * GraphDistance the length/total weight.
 *
 * Unreachable target: FindShortestPath -> {} (empty list), GraphDistance ->
 * Infinity.
 *
 * Memory (SPEC section 4): returns freshly-allocated results; frees res.
 */

#include "graph.h"
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include "assoc.h"   /* assoc_from_rules: the all-targets path table */
#include <stdlib.h>
#include <string.h>
#include <float.h>

/* BFS from src over out[]; fills parent[] (-1 = root/unvisited) and dist[]
 * (-1 = unreached). Caller allocates parent/dist of length a->n. */
static void bfs(const GraphAdj* a, int src, int* parent, int* dist) {
    for (int i = 0; i < a->n; i++) { parent[i] = -1; dist[i] = -1; }
    int* q = calloc((size_t)(a->n > 0 ? a->n : 1), sizeof(int));
    int head = 0, tail = 0;
    dist[src] = 0; q[tail++] = src;
    while (head < tail) {
        int u = q[head++];
        for (int j = 0; j < a->outdeg[u]; j++) {
            int w = a->out[u][j];
            if (dist[w] < 0) { dist[w] = dist[u] + 1; parent[w] = u; q[tail++] = w; }
        }
    }
    free(q);
}

/* ---- Weighted adjacency, local to this file (see header comment) --------- */

typedef struct {
    int n;
    const Expr* verts;   /* borrowed from g */
    int*    outdeg;
    int**   out;         /* out[v][k]   = neighbor vertex index */
    Expr*** outw;        /* outw[v][k]  = borrowed weight Expr* for that edge */
    Expr*   weights;     /* owns every Expr* referenced by outw; freed by wadj_free */
} WAdj;

static void wadj_free(WAdj* w) {
    if (!w) return;
    for (int i = 0; i < w->n; i++) { free(w->out[i]); free(w->outw[i]); }
    free(w->out); free(w->outw); free(w->outdeg);
    if (w->weights) expr_free(w->weights);
    free(w);
}

/* Builds a weighted successor adjacency directly from g's canonical form
 * (verts/edges/EdgeWeight), independent of GraphAdj. NULL if g is not a
 * valid, weights-usable graph. */
static WAdj* build_wadj(const Expr* g) {
    if (!graph_weights_usable(g)) return NULL;
    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    Expr* weights = graph_resolve_edge_weights(g);
    if (!weights) return NULL;

    WAdj* w = calloc(1, sizeof(WAdj));
    if (!w) { expr_free(weights); return NULL; }
    w->n = n; w->verts = verts; w->weights = weights;
    w->outdeg = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    w->out    = calloc((size_t)(n > 0 ? n : 1), sizeof(int*));
    w->outw   = calloc((size_t)(n > 0 ? n : 1), sizeof(Expr**));

    GraphVIdx* ix = graph_vidx_new((size_t)n);
    for (int i = 0; i < n; i++)
        graph_vidx_put(ix, verts->data.function.args[i], i);

    /* Pass 1: count out-degrees (each undirected edge contributes to both endpoints). */
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        const char* kind = graph_edge_kind(e);
        int ia = graph_vidx_get(ix, e->data.function.args[0]);
        int ib = graph_vidx_get(ix, e->data.function.args[1]);
        if (ia < 0 || ib < 0) continue;
        w->outdeg[ia]++;
        if (kind == SYM_UndirectedEdge) w->outdeg[ib]++;
    }
    for (int i = 0; i < n; i++) {
        w->out[i]  = w->outdeg[i] ? calloc((size_t)w->outdeg[i], sizeof(int)) : NULL;
        w->outw[i] = w->outdeg[i] ? calloc((size_t)w->outdeg[i], sizeof(Expr*)) : NULL;
    }
    int* fill = calloc((size_t)(n > 0 ? n : 1), sizeof(int));

    /* Pass 2: fill. */
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        const char* kind = graph_edge_kind(e);
        int ia = graph_vidx_get(ix, e->data.function.args[0]);
        int ib = graph_vidx_get(ix, e->data.function.args[1]);
        if (ia < 0 || ib < 0) continue;
        Expr* wt = weights->data.function.args[k];
        w->out[ia][fill[ia]]  = ib;
        w->outw[ia][fill[ia]] = wt;
        fill[ia]++;
        if (kind == SYM_UndirectedEdge) {
            w->out[ib][fill[ib]]  = ia;
            w->outw[ib][fill[ib]] = wt;
            fill[ib]++;
        }
    }
    free(fill);
    graph_vidx_free(ix);
    return w;
}

/* Dijkstra from src over w; fills parent[] (-1 = root/unvisited) and dist[]
 * (DBL_MAX = unreached). O(V^2) array scan, no heap -- consistent with this
 * subsystem's existing small-graph exact-algorithm precedent
 * (VertexConnectivity's own brute-force). dist[t] is also GraphDistance's
 * answer: a machine real, as in the Wolfram Language. */
static void dijkstra(const WAdj* w, int src, int* parent, double* dist) {
    char* done = calloc((size_t)(w->n > 0 ? w->n : 1), sizeof(char));
    for (int i = 0; i < w->n; i++) { parent[i] = -1; dist[i] = DBL_MAX; }
    dist[src] = 0.0;

    for (int iter = 0; iter < w->n; iter++) {
        int u = -1;
        double best = DBL_MAX;
        for (int i = 0; i < w->n; i++)
            if (!done[i] && dist[i] < best) { best = dist[i]; u = i; }
        if (u < 0) break;                 /* remaining vertices are unreachable */
        done[u] = 1;
        for (int j = 0; j < w->outdeg[u]; j++) {
            int v = w->out[u][j];
            double d = dist[u] + graph_weight_to_double(w->outw[u][j]);
            if (d < dist[v]) { dist[v] = d; parent[v] = u; }
        }
    }
    free(done);
}

/* Resolve g, s, t to a GraphAdj and endpoint indices. Returns adj (caller frees)
 * or NULL; on success sets *is,*it. */
/* True iff e is Method -> m with m Automatic or a shortest-path method name
 * the reference accepts.  Every method gives the same distances, so the
 * choice only selects an algorithm; the one here is Dijkstra (or BFS). */
static bool is_method_option(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2 ||
        e->data.function.head->type != EXPR_SYMBOL ||
        (e->data.function.head->data.symbol.name != SYM_Rule &&
         e->data.function.head->data.symbol.name != SYM_RuleDelayed)) return false;
    const Expr* k = e->data.function.args[0];
    const Expr* v = e->data.function.args[1];
    if (k->type != EXPR_SYMBOL || k->data.symbol.name != SYM_Method) return false;
    if (v->type == EXPR_SYMBOL) return v->data.symbol.name == SYM_Automatic;
    if (v->type != EXPR_STRING) return false;
    static const char* const M[] = { "Dijkstra", "BellmanFord", "FloydWarshall", "Johnson" };
    for (size_t i = 0; i < sizeof(M) / sizeof(M[0]); i++)
        if (strcmp(v->data.string, M[i]) == 0) return true;
    return false;
}

/* True when g carries an EdgeWeight list that the weighted search cannot use
 * (a symbolic or negative weight). */
static bool weights_unusable(const Expr* g) {
    const Expr* wl = graph_edge_weight_list(g);
    if (!wl || graph_weights_usable(g)) return false;
    /* only a symbolic weight: a negative number keeps the documented
     * hop-count fallback */
    for (size_t k = 0; k < wl->data.function.arg_count; k++)
        if (!expr_is_numeric_like(wl->data.function.args[k])) return true;
    return false;
}

static GraphAdj* resolve(Expr* res, int* is, int* it) {
    size_t argc = res->data.function.arg_count;
    if (argc == 4 && is_method_option(res->data.function.args[3])) argc = 3;
    if (argc != 3) return NULL;
    const Expr* g = res->data.function.args[0];
    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    *is = graph_vertex_index(a->verts, res->data.function.args[1]);
    *it = graph_vertex_index(a->verts, res->data.function.args[2]);
    if (*is < 0 || *it < 0) { graph_adj_free(a); return NULL; }
    return a;
}

/* ---- FindShortestPath[g, s, All] and FindShortestPath[g, All, All] ------ */

/* The shortest path from vertex `is` to every vertex of g, as the Association
 * <|v -> path, ...|> over VertexList order ({} for an unreachable v).  One
 * search (BFS, or Dijkstra under a usable EdgeWeight) serves every target, so
 * the all-targets form costs one search, not one per target.  NULL on
 * allocation failure. */
static Expr* sp_paths_from(const Expr* g, const GraphAdj* a, int is) {
    int n = a->n;
    int* parent = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    char* reach = calloc((size_t)(n > 0 ? n : 1), sizeof(char));
    Expr** rules = malloc(sizeof(Expr*) * (size_t)(n > 0 ? n : 1));
    int* path = malloc(sizeof(int) * (size_t)(n > 0 ? n : 1));
    bool ok = parent && reach && rules && path;
    if (ok && graph_weights_usable(g)) {
        WAdj* w = build_wadj(g);
        double* dist = w ? calloc((size_t)(n > 0 ? n : 1), sizeof(double)) : NULL;
        if (w && dist) {
            dijkstra(w, is, parent, dist);
            for (int v = 0; v < n; v++) reach[v] = dist[v] != DBL_MAX;
        } else ok = false;
        free(dist); wadj_free(w);
    } else if (ok) {
        int* dist = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
        if (dist) {
            bfs(a, is, parent, dist);
            for (int v = 0; v < n; v++) reach[v] = dist[v] >= 0;
        } else ok = false;
        free(dist);
    }
    Expr* out = NULL;
    if (ok) {
        for (int t = 0; t < n; t++) {
            int len = 0;
            if (reach[t]) {
                for (int v = t; ; v = parent[v]) { path[len++] = v; if (v == is) break; }
            }
            Expr** pv = malloc(sizeof(Expr*) * (size_t)(len > 0 ? len : 1));
            for (int k = 0; k < len; k++)
                pv[k] = expr_copy(a->verts->data.function.args[path[len - 1 - k]]);
            Expr* pl = expr_new_function(expr_new_symbol(SYM_List), pv, (size_t)len);
            free(pv);
            Expr* ra[2] = { expr_copy(a->verts->data.function.args[t]), pl };
            rules[t] = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
        }
        out = assoc_from_rules(rules, (size_t)n);
        for (int t = 0; t < n; t++) expr_free(rules[t]);
    }
    free(parent); free(reach); free(rules); free(path);
    return out;
}

/* Lookup[assoc, #k] (k = 1 or 2), adopting assoc. */
static Expr* sp_lookup_slot(Expr* assoc, int64_t k) {
    Expr* sa[1] = { expr_new_integer(k) };
    Expr* la[2] = { assoc, expr_new_function(expr_new_symbol(SYM_Slot), sa, 1) };
    return expr_new_function(expr_new_symbol(SYM_Lookup), la, 2);
}

/* FindShortestPath[g, s, All] is the function t |-> FindShortestPath[g, s, t]
 * and FindShortestPath[g, All, All] the function (s, t) |-> ... (audit
 * FindShortestPath#2).  Each is a pure function over a precomputed table of
 * paths, so applying it is a lookup: Function[Lookup[<|t -> path|>, #]] and
 * Function[Lookup[Lookup[<|s -> <|t -> path|>|>, #1], #2]].  NULL when the
 * call is not one of these forms or a source is not a vertex. */
static Expr* sp_all_targets(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc == 4 && is_method_option(res->data.function.args[3])) argc = 3;
    if (argc != 3) return NULL;
    const Expr* sv = res->data.function.args[1];
    const Expr* tv = res->data.function.args[2];
    if (!(tv->type == EXPR_SYMBOL && tv->data.symbol.name == SYM_All)) return NULL;
    bool all_src = sv->type == EXPR_SYMBOL && sv->data.symbol.name == SYM_All;
    const Expr* g = res->data.function.args[0];
    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    Expr* body = NULL;
    if (!all_src) {
        int is = graph_vertex_index(a->verts, sv);
        Expr* tab = is >= 0 ? sp_paths_from(g, a, is) : NULL;
        if (tab) body = sp_lookup_slot(tab, 1);
    } else {
        int n = a->n;
        Expr** rules = malloc(sizeof(Expr*) * (size_t)(n > 0 ? n : 1));
        int built = 0;
        for (int s0 = 0; rules && s0 < n; s0++) {
            Expr* tab = sp_paths_from(g, a, s0);
            if (!tab) break;
            Expr* ra[2] = { expr_copy(a->verts->data.function.args[s0]), tab };
            rules[built++] = expr_new_function(expr_new_symbol(SYM_Rule), ra, 2);
        }
        if (rules && built == n) {
            Expr* outer = assoc_from_rules(rules, (size_t)n);
            if (outer) {
                Expr* inner = sp_lookup_slot(outer, 1);
                Expr* sa[1] = { expr_new_integer(2) };
                Expr* la[2] = { inner, expr_new_function(expr_new_symbol(SYM_Slot), sa, 1) };
                body = expr_new_function(expr_new_symbol(SYM_Lookup), la, 2);
            }
        }
        for (int k = 0; k < built; k++) expr_free(rules[k]);
        free(rules);
    }
    graph_adj_free(a);
    if (!body) return NULL;
    return expr_new_function(expr_new_symbol(SYM_Function), &body, 1);
}

Expr* builtin_find_shortest_path(Expr* res) {
    {
        Expr* all = sp_all_targets(res);
        if (all) return all;
    }
    int is, it;
    GraphAdj* a = resolve(res, &is, &it);
    if (!a) return NULL;
    const Expr* g = res->data.function.args[0];

    Expr* out;
    if (graph_weights_usable(g)) {
        WAdj* w = build_wadj(g);
        int* parent = w ? calloc((size_t)(w->n > 0 ? w->n : 1), sizeof(int)) : NULL;
        double* dist = w ? calloc((size_t)(w->n > 0 ? w->n : 1), sizeof(double)) : NULL;
        if (!parent || !dist) {                 /* allocation failure: unevaluated */
            free(parent); free(dist); wadj_free(w); graph_adj_free(a);
            return NULL;
        }
        dijkstra(w, is, parent, dist);
        if (dist[it] == DBL_MAX) {
            out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        } else {
            /* Reconstruct the vertex path from parent[] (same shape as BFS's). */
            int len = 1, v = it;
            while (v != is) { len++; v = parent[v]; }
            Expr** path = calloc((size_t)len, sizeof(Expr*));
            v = it;
            for (int k = len - 1; k >= 0; k--) { path[k] = expr_copy(w->verts->data.function.args[v]); v = parent[v]; }
            out = expr_new_function(expr_new_symbol(SYM_List), path, (size_t)len);
            free(path);
        }
        free(parent); free(dist); wadj_free(w);
    } else {
        int* parent = calloc((size_t)a->n, sizeof(int));
        int* dist   = calloc((size_t)a->n, sizeof(int));
        bfs(a, is, parent, dist);
        if (dist[it] < 0) {
            out = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);  /* {} */
        } else {
            int len = dist[it] + 1;
            Expr** path = calloc((size_t)len, sizeof(Expr*));
            int v = it;
            for (int k = len - 1; k >= 0; k--) { path[k] = expr_copy(a->verts->data.function.args[v]); v = parent[v]; }
            out = expr_new_function(expr_new_symbol(SYM_List), path, (size_t)len);
            free(path);
        }
        free(parent); free(dist);
    }
    graph_adj_free(a);
    return out;
}

Expr* builtin_graph_distance(Expr* res) {
    /* A distance under a symbolic weight is not a hop count:
     * unevaluated, as GraphDistance[g, s] is.  (FindShortestPath keeps the
     * unweighted path: any path it gives is a path of g.) */
    if (res->data.function.arg_count >= 1 && weights_unusable(res->data.function.args[0]))
        return NULL;
    int is, it;
    GraphAdj* a = resolve(res, &is, &it);
    if (!a) return NULL;
    const Expr* g = res->data.function.args[0];

    Expr* out;
    if (graph_weights_usable(g)) {
        WAdj* w = build_wadj(g);
        int* parent = w ? calloc((size_t)(w->n > 0 ? w->n : 1), sizeof(int)) : NULL;
        double* dist = w ? calloc((size_t)(w->n > 0 ? w->n : 1), sizeof(double)) : NULL;
        if (!parent || !dist) {                 /* allocation failure: unevaluated */
            free(parent); free(dist); wadj_free(w); graph_adj_free(a);
            return NULL;
        }
        dijkstra(w, is, parent, dist);
        /* A weighted distance is a machine real, as in the Wolfram Language
         * (GraphDistance[g, 1, 3] with weights {5, 7} is 12., not 12): the
         * Dijkstra accumulator itself, summed along the path in path order, so
         * it agrees bit-for-bit with GraphDistance[g, s] / GraphDistanceMatrix. */
        out = (dist[it] == DBL_MAX) ? expr_new_symbol(SYM_Infinity)
                                    : expr_new_real(dist[it]);
        free(parent); free(dist); wadj_free(w);
    } else {
        int* parent = calloc((size_t)a->n, sizeof(int));
        int* dist   = calloc((size_t)a->n, sizeof(int));
        bfs(a, is, parent, dist);
        out = (dist[it] < 0) ? expr_new_symbol(SYM_Infinity)
                             : expr_new_integer(dist[it]);
        free(parent); free(dist);
    }
    graph_adj_free(a);
    return out;
}
