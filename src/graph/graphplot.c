/* graphplot.c - GraphPlot[g] and the shared graph renderer.
 *
 * graph_render() turns a validated graph into a Graphics[...] tree of the same
 * primitives the plotting engine already draws (Line, Disk, Text), preceded by
 * RGBColor directives for per-element styling. The notebook Plotly serializer
 * and the Raylib renderer both consume this with no special-casing, and the
 * text placeholder is used when USE_GRAPHICS=0.
 *
 * Vertex positions come from src/graph/layout.c, selected by the GraphLayout
 * option (circular by default). GraphPlot[g, opts] parses GraphLayout,
 * VertexStyle, EdgeStyle, VertexLabels and VertexSize; HighlightGraph reuses
 * graph_render() through the highlight masks in GraphStyle.
 *
 * Memory (SPEC section 4): builtins return a freshly-allocated Graphics tree;
 * the evaluator frees res. graph_render is read-only over its arguments.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include "eval.h"
#include <math.h>
#include <stdlib.h>

#define NODE_RADIUS 0.08

/* Default palette (RGB in 0..1). */
static const double DEF_VERTEX[3] = { 0.20, 0.42, 0.82 };  /* blue             */
static const double DEF_EDGE[3]   = { 0.52, 0.56, 0.66 };  /* slate            */
static const double ACCENT[3]     = { 0.95, 0.45, 0.10 };  /* highlight orange */
static const double DIM[3]        = { 0.80, 0.83, 0.88 };  /* dimmed grey      */

static Expr* point2(double x, double y) {
    Expr* xy[2] = { expr_new_real(x), expr_new_real(y) };
    return expr_new_function(expr_new_symbol(SYM_List), xy, 2);
}

static Expr* rgb(const double c[3]) {
    Expr* a[3] = { expr_new_real(c[0]), expr_new_real(c[1]), expr_new_real(c[2]) };
    return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
}

/* A color directive for element `idx`, honoring highlight masks first, then the
 * explicit style color, then the default. `mask` may be NULL. */
static Expr* color_for(const char* mask, int idx, const Expr* style_color,
                       const double def[3]) {
    if (mask) return rgb(mask[idx] ? ACCENT : DIM);
    if (style_color) return expr_copy((Expr*)style_color);
    return rgb(def);
}

Expr* graph_render(const Expr* g, const GraphStyle* st) {
    if (!graph_is_valid(g)) return NULL;
    GraphStyle defaults = {0};
    defaults.show_labels = 1;
    if (!st) st = &defaults;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;
    double radius = (st->vertex_size > 0.0) ? st->vertex_size : NODE_RADIUS;
    int labels = st->show_labels;

    double* x = (n > 0) ? calloc((size_t)n, sizeof(double)) : NULL;
    double* y = (n > 0) ? calloc((size_t)n, sizeof(double)) : NULL;
    if (n > 0 && (!x || !y)) { free(x); free(y); return NULL; }
    graph_compute_layout(g, st->layout, x, y);

    /* Worst case: (color + line) per edge, (color + disk) per vertex, and one
     * label per vertex. */
    size_t cap = ne * 2 + (size_t)n * 3 + 1;
    Expr** prims = (cap > 0) ? calloc(cap, sizeof(Expr*)) : NULL;
    size_t p = 0;

    /* Edges first so vertices draw on top. */
    for (size_t k = 0; k < ne; k++) {
        const Expr* e = edges->data.function.args[k];
        int iu = graph_vertex_index(verts, e->data.function.args[0]);
        int iv = graph_vertex_index(verts, e->data.function.args[1]);
        if (iu < 0 || iv < 0) continue;
        prims[p++] = color_for(st->hi_edge, (int)k, st->edge_color, DEF_EDGE);
        Expr* pts[2] = { point2(x[iu], y[iu]), point2(x[iv], y[iv]) };
        Expr* seg = expr_new_function(expr_new_symbol(SYM_List), pts, 2);
        Expr* la[1] = { seg };
        prims[p++] = expr_new_function(expr_new_symbol(SYM_Line), la, 1);
    }
    for (int i = 0; i < n; i++) {
        prims[p++] = color_for(st->hi_vert, i, st->vertex_color, DEF_VERTEX);
        Expr* da[2] = { point2(x[i], y[i]), expr_new_real(radius) };
        prims[p++] = expr_new_function(expr_new_symbol(SYM_Disk), da, 2);
    }
    if (labels) {
        for (int i = 0; i < n; i++) {
            Expr* ta[2] = { expr_copy((Expr*)verts->data.function.args[i]), point2(x[i], y[i]) };
            prims[p++] = expr_new_function(expr_new_symbol(SYM_Text), ta, 2);
        }
    }

    free(x); free(y);
    Expr* prim_list = expr_new_function(expr_new_symbol(SYM_List), prims, p);
    free(prims);
    Expr* ga[1] = { prim_list };
    return expr_new_function(expr_new_symbol(SYM_Graphics), ga, 1);
}

Expr* graph_default_graphics(const Expr* g) {
    return graph_render(g, NULL);
}

/* ---- option parsing ------------------------------------------------------- */

/* Resolve a color option value to an owned RGBColor[...] (converting GrayLevel
 * and named colors, which evaluate to RGBColor). Returns NULL if not a color. */
static Expr* resolve_color(const Expr* v) {
    Expr* ev = evaluate(expr_copy((Expr*)v));
    if (!ev) return NULL;
    if (ev->type == EXPR_FUNCTION && ev->data.function.head
        && ev->data.function.head->type == EXPR_SYMBOL) {
        const char* h = ev->data.function.head->data.symbol;
        if (h == SYM_RGBColor && ev->data.function.arg_count >= 3) return ev;
        if (h == SYM_GrayLevel && ev->data.function.arg_count >= 1) {
            Expr* lv = ev->data.function.args[0];
            Expr* a[3] = { expr_copy(lv), expr_copy(lv), expr_copy(lv) };
            expr_free(ev);
            return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
        }
    }
    expr_free(ev);
    return NULL;
}

/* GraphPlot[g] / GraphPlot[g, opts...] */
Expr* builtin_graph_plot(Expr* res) {
    if (res->data.function.arg_count < 1) return NULL;
    const Expr* g = res->data.function.args[0];
    if (!graph_is_valid(g)) return NULL;

    GraphStyle st = {0};
    st.show_labels = 1;
    Expr* vcol = NULL; Expr* ecol = NULL;   /* owned; freed before return       */

    for (size_t i = 1; i < res->data.function.arg_count; i++) {
        const Expr* opt = res->data.function.args[i];
        if (opt->type != EXPR_FUNCTION || opt->data.function.arg_count != 2
            || !opt->data.function.head
            || opt->data.function.head->type != EXPR_SYMBOL
            || !(opt->data.function.head->data.symbol == SYM_Rule
                 || opt->data.function.head->data.symbol == SYM_RuleDelayed))
            continue;
        const char* name = (opt->data.function.args[0]->type == EXPR_SYMBOL)
            ? opt->data.function.args[0]->data.symbol : NULL;
        const Expr* rhs = opt->data.function.args[1];
        if (name == SYM_GraphLayout) {
            if (rhs->type == EXPR_STRING) st.layout = rhs->data.string;
        } else if (name == SYM_VertexStyle) {
            if (vcol) expr_free(vcol);
            vcol = resolve_color(rhs); st.vertex_color = vcol;
        } else if (name == SYM_EdgeStyle) {
            if (ecol) expr_free(ecol);
            ecol = resolve_color(rhs); st.edge_color = ecol;
        } else if (name == SYM_VertexSize) {
            Expr* ev = evaluate(expr_copy((Expr*)rhs));
            if (ev && ev->type == EXPR_REAL)         st.vertex_size = ev->data.real;
            else if (ev && ev->type == EXPR_INTEGER) st.vertex_size = (double)ev->data.integer;
            if (ev) expr_free(ev);
        } else if (name == SYM_VertexLabels) {
            if (rhs->type == EXPR_SYMBOL
                && (rhs->data.symbol == SYM_None || rhs->data.symbol == SYM_False))
                st.show_labels = 0;
        }
    }

    Expr* out = graph_render(g, &st);
    if (vcol) expr_free(vcol);
    if (ecol) expr_free(ecol);
    return out;
}
