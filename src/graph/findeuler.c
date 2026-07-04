/* findeuler.c - FindEulerianCycle[g]: an Eulerian cycle (a closed walk using
 * every edge exactly once) as a vertex list {v0, v1, ..., v0}, or {} if none
 * exists. The constructive companion to EulerianGraphQ.
 *
 * Hierholzer's algorithm, O(V+E): walk unused edges onto a stack, splicing in
 * sub-tours as dead ends are hit; the reversed pop-order is the Euler tour. It
 * doubles as the existence test — if the tour fails to use all ne edges (odd
 * degree / in!=out / disconnected), the graph is not Eulerian and we return {}.
 * Works for directed (follow out-edges) and undirected (each edge usable once
 * from either end).
 *
 * Memory (SPEC section 4): returns a fresh List; frees res. NULL on a non-graph.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

Expr* builtin_find_eulerian_cycle(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    int ne = (int)edges->data.function.arg_count;
    Expr* empty = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    if (ne == 0) return empty;   /* no edges -> no cycle */

    int directed = 1;
    for (int i = 0; i < ne; i++)
        if (graph_edge_kind(edges->data.function.args[i]) != SYM_DirectedEdge) { directed = 0; break; }

    /* Incidence lists with edge ids (each undirected edge appears from both
     * ends but shares one id, so it is walked at most once). */
    int* deg = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    int ia_b, ib_b;
    for (int k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        ia_b = graph_vertex_index(verts, e->data.function.args[0]);
        ib_b = graph_vertex_index(verts, e->data.function.args[1]);
        deg[ia_b]++; if (!directed) deg[ib_b]++;
    }
    int** nbr = calloc((size_t)(n > 0 ? n : 1), sizeof(int*));
    int** eid = calloc((size_t)(n > 0 ? n : 1), sizeof(int*));
    int* pos = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    for (int v = 0; v < n; v++) {
        nbr[v] = (deg[v] > 0) ? malloc((size_t)deg[v] * sizeof(int)) : NULL;
        eid[v] = (deg[v] > 0) ? malloc((size_t)deg[v] * sizeof(int)) : NULL;
    }
    for (int k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        int a = graph_vertex_index(verts, e->data.function.args[0]);
        int b = graph_vertex_index(verts, e->data.function.args[1]);
        nbr[a][pos[a]] = b; eid[a][pos[a]++] = k;
        if (!directed) { nbr[b][pos[b]] = a; eid[b][pos[b]++] = k; }
    }

    char* used = calloc((size_t)ne, 1);
    int* cursor = calloc((size_t)(n > 0 ? n : 1), sizeof(int));
    int start = 0;
    for (int v = 0; v < n; v++) if (deg[v] > 0) { start = v; break; }

    int* stack = malloc((size_t)(ne + 1) * sizeof(int));
    int* circuit = malloc((size_t)(ne + 1) * sizeof(int));
    int sp = 0, cp = 0;
    stack[sp++] = start;
    while (sp > 0) {
        int v = stack[sp - 1];
        while (cursor[v] < deg[v] && used[eid[v][cursor[v]]]) cursor[v]++;
        if (cursor[v] < deg[v]) {
            int idx = cursor[v]++;
            used[eid[v][idx]] = 1;
            stack[sp++] = nbr[v][idx];
        } else {
            circuit[cp++] = v;
            sp--;
        }
    }

    /* An Eulerian CYCLE uses every edge exactly once AND is closed. Requiring
     * both distinguishes it from an Eulerian path (all edges, open walk) and
     * from a traversal that stalled before covering the graph. */
    int nused = 0;
    for (int k = 0; k < ne; k++) nused += used[k];
    int closed = (cp >= 2 && circuit[0] == circuit[cp - 1]);

    Expr* out;
    if (nused == ne && cp == ne + 1 && closed) {
        Expr** items = calloc((size_t)cp, sizeof(Expr*));
        for (int i = 0; i < cp; i++)     /* circuit is reversed */
            items[i] = expr_copy((Expr*)verts->data.function.args[circuit[cp - 1 - i]]);
        out = expr_new_function(expr_new_symbol(SYM_List), items, (size_t)cp);
        free(items);
        expr_free(empty);
    } else {
        out = empty;                     /* not Eulerian */
    }

    for (int v = 0; v < n; v++) { free(nbr[v]); free(eid[v]); }
    free(nbr); free(eid); free(pos); free(deg); free(used); free(cursor);
    free(stack); free(circuit);
    return out;
}
