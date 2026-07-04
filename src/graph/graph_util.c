/* graph_util.c - shared, read-only helpers for the graph subsystem.
 *
 * Graphs are ordinary Expr trees; these helpers inspect the canonical form
 *
 *     Graph[ List[v1, ...], List[edge1, ...] ]
 *
 * where each edge is a 2-argument DirectedEdge[u, v] or UndirectedEdge[u, v].
 * Nothing here allocates or mutates; ownership contracts live in the callers.
 *
 * The MVP vertex-membership test is a linear expr_eq scan (O(V) per lookup).
 * That is fine for pico-CAS graph sizes; an expr_hash-based index is the
 * documented upgrade path when profiling warrants it.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* True iff e is a function node whose head is the interned symbol `sym`. */
static int head_is_sym(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == sym;
}

int graph_is_list(const Expr* e) {
    return head_is_sym(e, SYM_List);
}

const char* graph_edge_kind(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 2)
        return NULL;
    if (head_is_sym(e, SYM_DirectedEdge))   return SYM_DirectedEdge;
    if (head_is_sym(e, SYM_UndirectedEdge)) return SYM_UndirectedEdge;
    return NULL;
}

/* True iff `v` is structurally equal to some element of List `list`. */
static int vertex_in_list(const Expr* list, const Expr* v) {
    for (size_t i = 0; i < list->data.function.arg_count; i++) {
        if (expr_eq(list->data.function.args[i], v)) return 1;
    }
    return 0;
}

int graph_vertex_index(const Expr* verts, const Expr* v) {
    if (!graph_is_list(verts)) return -1;
    for (size_t i = 0; i < verts->data.function.arg_count; i++) {
        if (expr_eq(verts->data.function.args[i], v)) return (int)i;
    }
    return -1;
}

/* ---- Phase 5: adjacency scaffolding --------------------------------------- */

void graph_adj_free(GraphAdj* a) {
    if (!a) return;
    for (int i = 0; i < a->n; i++) { free(a->out[i]); free(a->in[i]); }
    free(a->out); free(a->in);
    free(a->outdeg); free(a->indeg);
    free(a);
}

GraphAdj* graph_build_adj(const Expr* g) {
    if (!graph_is_valid(g)) return NULL;
    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    GraphAdj* a = calloc(1, sizeof(GraphAdj));
    if (!a) return NULL;
    a->n = n;
    a->verts = verts;
    a->outdeg = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    a->indeg  = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    a->out    = calloc((size_t)(n > 0 ? n : 1), sizeof(int*));
    a->in     = calloc((size_t)(n > 0 ? n : 1), sizeof(int*));

    /* Pass 1: count degrees. */
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        const char* kind = graph_edge_kind(e);
        int ia = graph_vertex_index(verts, e->data.function.args[0]);
        int ib = graph_vertex_index(verts, e->data.function.args[1]);
        a->outdeg[ia]++; a->indeg[ib]++;
        if (kind == SYM_UndirectedEdge) { a->outdeg[ib]++; a->indeg[ia]++; }
    }
    for (int i = 0; i < n; i++) {
        a->out[i] = (a->outdeg[i] > 0) ? calloc((size_t)a->outdeg[i], sizeof(int)) : NULL;
        a->in[i]  = (a->indeg[i]  > 0) ? calloc((size_t)a->indeg[i],  sizeof(int)) : NULL;
    }

    /* Pass 2: fill (reuse degree counters as write cursors). */
    int* oc = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    int* ic = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        const char* kind = graph_edge_kind(e);
        int ia = graph_vertex_index(verts, e->data.function.args[0]);
        int ib = graph_vertex_index(verts, e->data.function.args[1]);
        a->out[ia][oc[ia]++] = ib;  a->in[ib][ic[ib]++] = ia;
        if (kind == SYM_UndirectedEdge) {
            a->out[ib][oc[ib]++] = ia;  a->in[ia][ic[ia]++] = ib;
        }
    }
    free(oc); free(ic);
    return a;
}

int graph_count_components(const GraphAdj* a, const char* removed, int* active_out) {
    int n = a->n;
    char* seen = calloc((size_t)(n > 0 ? n : 1), sizeof(char));
    int* stack = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    int comps = 0, active = 0;

    for (int s = 0; s < n; s++) {
        if (removed && removed[s]) continue;
        active++;
        if (seen[s]) continue;
        /* New component: DFS over underlying undirected neighbors (out + in). */
        comps++;
        int top = 0; stack[top++] = s; seen[s] = 1;
        while (top > 0) {
            int u = stack[--top];
            for (int j = 0; j < a->outdeg[u]; j++) {
                int w = a->out[u][j];
                if ((removed && removed[w]) || seen[w]) continue;
                seen[w] = 1; stack[top++] = w;
            }
            for (int j = 0; j < a->indeg[u]; j++) {
                int w = a->in[u][j];
                if ((removed && removed[w]) || seen[w]) continue;
                seen[w] = 1; stack[top++] = w;
            }
        }
    }
    free(seen); free(stack);
    if (active_out) *active_out = active;
    return comps;
}

/* Two normalized edges are "parallel" (duplicates) when they connect the same
 * endpoints in a way the graph cannot distinguish:
 *   - directed:   same head and same ordered (u, v);
 *   - undirected: same head and the same unordered pair {u, v}.
 * Directed a->b and b->a are distinct; that is allowed. */
static int edges_parallel(const Expr* e1, const Expr* e2) {
    const char* k1 = graph_edge_kind(e1);
    const char* k2 = graph_edge_kind(e2);
    if (!k1 || k1 != k2) return 0;
    const Expr* a1 = e1->data.function.args[0];
    const Expr* b1 = e1->data.function.args[1];
    const Expr* a2 = e2->data.function.args[0];
    const Expr* b2 = e2->data.function.args[1];
    if (expr_eq(a1, a2) && expr_eq(b1, b2)) return 1;
    if (k1 == SYM_UndirectedEdge && expr_eq(a1, b2) && expr_eq(b1, a2)) return 1;
    return 0;
}

int graph_is_valid_head(const Expr* g, const char* head_sym) {
    if (!head_is_sym(g, head_sym) || g->data.function.arg_count != 2)
        return 0;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    if (!graph_is_list(verts) || !graph_is_list(edges)) return 0;

    size_t ne = edges->data.function.arg_count;
    for (size_t i = 0; i < ne; i++) {
        const Expr* edge = edges->data.function.args[i];
        if (!graph_edge_kind(edge)) return 0;              /* un-normalized/3-arg */

        const Expr* u = edge->data.function.args[0];
        const Expr* v = edge->data.function.args[1];
        if (expr_eq(u, v)) return 0;                        /* self-loop */
        if (!vertex_in_list(verts, u)) return 0;            /* endpoint not a vertex */
        if (!vertex_in_list(verts, v)) return 0;

        for (size_t j = 0; j < i; j++) {
            if (edges_parallel(edge, edges->data.function.args[j])) return 0;
        }
    }
    return 1;
}

int graph_is_valid(const Expr* g) {
    /* Both 2D Graph and 3D Graph3D count as valid graphs, so every query,
     * matrix, and algorithm builtin works transparently on a Graph3D. */
    return graph_is_valid_head(g, SYM_Graph)
        || graph_is_valid_head(g, SYM_Graph3D);
}
