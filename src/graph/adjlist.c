/* adjlist.c - AdjacencyList[g] and AdjacencyList[g, v].
 *
 * Neighbors of v: successors for directed edges (v -> u yields u), and both
 * endpoints for undirected edges. This matches the row convention of
 * AdjacencyMatrix (a 1 in row v, column u means v is adjacent to u). Neighbors
 * are returned in first-appearance order, de-duplicated.
 *
 * AdjacencyList[g]     -> {neighbors(v1), neighbors(v2), ...} in vertex order.
 * AdjacencyList[g, v]  -> neighbors(v).
 * AdjacencyList[g, patt] -> the neighbours of every vertex matching patt,
 *                         in VertexList order.
 *
 * Memory (SPEC section 4): returns freshly-allocated lists; evaluator frees res.
 */

#include "graph.h"
#include "graph_ops.h"   /* gops_has_pattern, gops_matchq */
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <stdint.h>

/* Build List of neighbors of `v` (deduplicated, first-appearance order). */
static Expr* neighbors_of(const Expr* g, const Expr* v) {
    const Expr* edges = g->data.function.args[1];
    size_t ne = edges->data.function.arg_count;
    Expr** nbr = (ne > 0) ? calloc(ne * 2, sizeof(Expr*)) : NULL;
    size_t n = 0;

    for (size_t i = 0; i < ne; i++) {
        const Expr* e = edges->data.function.args[i];
        const char* kind = graph_edge_kind(e);
        const Expr* a = e->data.function.args[0];
        const Expr* b = e->data.function.args[1];
        const Expr* add = NULL;
        if (kind == SYM_UndirectedEdge) {
            if (expr_eq(a, v)) add = b;
            else if (expr_eq(b, v)) add = a;
        } else {
            if (expr_eq(a, v)) add = b;   /* successor only */
        }
        if (add) {
            int seen = 0;
            for (size_t j = 0; j < n; j++)
                if (expr_eq(nbr[j], add)) { seen = 1; break; }
            if (!seen) nbr[n++] = expr_copy((Expr*)add);
        }
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), nbr, n);
    free(nbr);
    return list;
}

/* AdjacencyList[g, patt]: the vertices adjacent to some vertex matching
 * patt, each once, in VertexList order. */
static Expr* neighbors_of_pattern(const Expr* g, const Expr* patt) {
    const Expr* verts = g->data.function.args[0];
    size_t nv = verts->data.function.arg_count;
    unsigned char* hit = calloc(nv ? nv : 1, 1);
    if (!hit) return NULL;
    for (size_t i = 0; i < nv; i++) {
        Expr* vi = verts->data.function.args[i];
        if (!gops_matchq(vi, (Expr*)patt)) continue;
        Expr* nb = neighbors_of(g, vi);
        for (size_t j = 0; nb && j < nb->data.function.arg_count; j++) {
            int q = graph_vertex_position(g, nb->data.function.args[j]);
            if (q >= 0) hit[q] = 1;
        }
        expr_free(nb);
    }
    size_t m = 0;
    Expr** items = malloc(sizeof(Expr*) * (nv ? nv : 1));
    for (size_t i = 0; items && i < nv; i++)
        if (hit[i]) items[m++] = expr_copy(verts->data.function.args[i]);
    free(hit);
    if (!items) return NULL;
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, m);
    free(items);
    return out;
}

/* AdjacencyList[g, v, d]: the vertices other than v within d steps of v
 * (following the same out-neighbourhood as the two-argument form), in
 * VertexList order. d is a nonnegative Integer or Infinity. */
static Expr* within_distance(const Expr* g, int src, int64_t d) {
    GraphAdj* a = graph_build_adj(g);
    if (!a) return NULL;
    int n = a->n;
    int64_t* dist = malloc(sizeof(int64_t) * (size_t)(n > 0 ? n : 1));
    int* queue = malloc(sizeof(int) * (size_t)(n > 0 ? n : 1));
    Expr* out = NULL;
    if (dist && queue) {
        for (int i = 0; i < n; i++) dist[i] = -1;
        int head = 0, tail = 0;
        dist[src] = 0; queue[tail++] = src;
        while (head < tail) {
            int u = queue[head++];
            if (dist[u] >= d) continue;
            for (int k = 0; k < a->outdeg[u]; k++) {
                int w = a->out[u][k];
                if (dist[w] < 0) { dist[w] = dist[u] + 1; queue[tail++] = w; }
            }
        }
        Expr** items = malloc(sizeof(Expr*) * (size_t)(n > 0 ? n : 1));
        size_t m = 0;
        for (int i = 0; items && i < n; i++)
            if (i != src && dist[i] >= 0)
                items[m++] = expr_copy(a->verts->data.function.args[i]);
        if (items) {
            out = expr_new_function(expr_new_symbol(SYM_List), items, m);
            free(items);
        }
    }
    free(dist); free(queue);
    graph_adj_free(a);
    return out;
}

Expr* builtin_adjacency_list(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 3) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;
    const Expr* verts = g->data.function.args[0];

    if (argc == 3) {
        const Expr* de = res->data.function.args[2];
        int64_t d;
        if (de->type == EXPR_INTEGER && de->data.integer >= 0) d = de->data.integer;
        else if (de->type == EXPR_SYMBOL && de->data.symbol.name == SYM_Infinity) d = INT64_MAX;
        else return NULL;
        int src = graph_vertex_position(g, res->data.function.args[1]);
        if (src < 0) return NULL;
        return within_distance(g, src, d);
    }

    if (argc == 2) {
        const Expr* v = res->data.function.args[1];
        if (graph_vertex_index(verts, v) >= 0) return neighbors_of(g, v);
        if (!gops_has_pattern(v)) return NULL;
        return neighbors_of_pattern(g, v);
    }

    size_t nv = verts->data.function.arg_count;
    Expr** rows = (nv > 0) ? calloc(nv, sizeof(Expr*)) : NULL;
    if (nv > 0 && !rows) return NULL;
    for (size_t i = 0; i < nv; i++)
        rows[i] = neighbors_of(g, verts->data.function.args[i]);
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), rows, nv);
    free(rows);
    return list;
}
