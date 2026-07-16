/* graphplot.c - GraphPlot[g]: render a graph as a Graphics[...] expression.
 *
 * Emits the same primitives the plotting engine uses (Line, Disk, Text), so the
 * existing renderer draws it with no renderer changes (and the text placeholder
 * is used when USE_GRAPHICS=0). Vertices are laid out on a circle (MVP layout;
 * a force-directed spring layout is the documented future hook). Each edge is a
 * Line between its endpoints, each vertex a Disk plus a Text label.
 *
 * Directed edges are drawn as plain lines in the MVP (no arrowheads yet).
 *
 * Memory (SPEC section 4): returns a freshly-allocated Graphics tree; the
 * evaluator frees res.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define NODE_RADIUS 0.08

static Expr* point2(double x, double y) {
    Expr* xy[2] = { expr_new_real(x), expr_new_real(y) };
    return expr_new_function(expr_new_symbol(SYM_List), xy, 2);
}

Expr* builtin_graph_plot(Expr* res) {
    if (res->data.function.arg_count != 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    /* Circular layout coordinates. */
    double* x = (n > 0) ? calloc((size_t)n, sizeof(double)) : NULL;
    double* y = (n > 0) ? calloc((size_t)n, sizeof(double)) : NULL;
    for (int i = 0; i < n; i++) {
        if (n == 1) { x[i] = 0.0; y[i] = 0.0; }
        else {
            double t = 2.0 * M_PI * (double)i / (double)n;
            x[i] = cos(t); y[i] = sin(t);
        }
    }

    /* Primitives: one Line per edge, then one Disk + one Text per vertex. */
    size_t total = ne + (size_t)n * 2;
    Expr** prims = (total > 0) ? calloc(total, sizeof(Expr*)) : NULL;
    size_t p = 0;

    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        int iu = graph_vertex_index(verts, e->data.function.args[0]);
        int iv = graph_vertex_index(verts, e->data.function.args[1]);
        Expr* pts[2] = { point2(x[iu], y[iu]), point2(x[iv], y[iv]) };
        Expr* seg = expr_new_function(expr_new_symbol(SYM_List), pts, 2);
        Expr* la[1] = { seg };
        prims[p++] = expr_new_function(expr_new_symbol(SYM_Line), la, 1);
    }
    for (int i = 0; i < n; i++) {
        Expr* da[2] = { point2(x[i], y[i]), expr_new_real(NODE_RADIUS) };
        prims[p++] = expr_new_function(expr_new_symbol(SYM_Disk), da, 2);
    }
    for (int i = 0; i < n; i++) {
        Expr* ta[2] = { expr_copy(verts->data.function.args[i]), point2(x[i], y[i]) };
        prims[p++] = expr_new_function(expr_new_symbol(SYM_Text), ta, 2);
    }

    free(x); free(y);
    Expr* prim_list = expr_new_function(expr_new_symbol(SYM_List), prims, total);
    free(prims);
    Expr* ga[1] = { prim_list };
    return expr_new_function(expr_new_symbol(SYM_Graphics), ga, 1);
}
