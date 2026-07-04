/* linegraph.c - LineGraph[g]: the line graph L(g).
 *
 * The vertices of L(g) are the edges of g; two are adjacent when the underlying
 * edges meet:
 *   undirected g -> undirected L(g): edges {a,b} and {c,d} are joined iff they
 *                   share an endpoint.
 *   directed g   -> directed L(g): an arc a->b points to an arc c->d iff b == c
 *                   (head of one is the tail of the next).
 * Line-graph vertices are the edge expressions themselves (e.g.
 * UndirectedEdge[1,2]); a simple graph has distinct edges, so they form a valid
 * vertex set.
 *
 * O(E^2) pairwise endpoint comparison. Memory (SPEC section 4): returns the
 * canonical Graph directly; frees res. NULL on a non-graph argument.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

/* True if undirected edges e,f share an endpoint. */
static int shares_endpoint(const Expr* e, const Expr* f) {
    Expr* a = e->data.function.args[0]; Expr* b = e->data.function.args[1];
    Expr* c = f->data.function.args[0]; Expr* d = f->data.function.args[1];
    return expr_eq(a, c) || expr_eq(a, d) || expr_eq(b, c) || expr_eq(b, d);
}

Expr* builtin_line_graph(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* edges = g->data.function.args[1];
    size_t ne = edges->data.function.arg_count;

    int directed = (ne > 0);
    for (size_t i = 0; i < ne; i++)
        if (graph_edge_kind(edges->data.function.args[i]) != SYM_DirectedEdge) { directed = 0; break; }

    /* Line-graph vertices = copies of the edges. */
    Expr** lverts = (ne > 0) ? calloc(ne, sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < ne; i++) lverts[i] = expr_copy(edges->data.function.args[i]);

    size_t cap = directed ? ne * (ne > 0 ? ne - 1 : 0) : ne * (ne > 0 ? ne - 1 : 0) / 2;
    Expr** ledges = (cap > 0) ? calloc(cap, sizeof(Expr*)) : NULL;
    size_t m = 0;

    if (directed) {
        for (size_t i = 0; i < ne; i++) {
            const Expr* ei = edges->data.function.args[i];
            for (size_t j = 0; j < ne; j++) {
                if (i == j) continue;
                const Expr* ej = edges->data.function.args[j];
                if (expr_eq(ei->data.function.args[1], ej->data.function.args[0])) {
                    Expr* a[2] = { expr_copy((Expr*)ei), expr_copy((Expr*)ej) };
                    ledges[m++] = expr_new_function(expr_new_symbol(SYM_DirectedEdge), a, 2);
                }
            }
        }
    } else {
        for (size_t i = 0; i < ne; i++)
            for (size_t j = i + 1; j < ne; j++)
                if (shares_endpoint(edges->data.function.args[i], edges->data.function.args[j])) {
                    Expr* a[2] = { expr_copy(edges->data.function.args[i]),
                                   expr_copy(edges->data.function.args[j]) };
                    ledges[m++] = expr_new_function(expr_new_symbol(SYM_UndirectedEdge), a, 2);
                }
    }

    Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), lverts, ne);
    Expr* elist = expr_new_function(expr_new_symbol(SYM_List), ledges, m);
    free(lverts); free(ledges);
    Expr* gargs[2] = { vlist, elist };
    return expr_new_function(expr_new_symbol(SYM_Graph), gargs, 2);
}
