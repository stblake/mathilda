/* render3d.c - render a Graph3D as a Graphics3D[...] node-link diagram.
 *
 * Mirrors graph_render() (src/graph/graphplot.c) but in three dimensions: edges
 * become 3D Lines and vertices a Point set, each preceded by an RGBColor
 * directive. The notebook serializer (src/graphics/graphics_json.c) turns a
 * Graphics3D into Plotly scatter3d traces; the REPL auto-displays a bare
 * Graph3D the same way it does Graph / Graphics.
 *
 * Coordinates come from graph_compute_layout3d() (src/graph/layout.c). Vertex
 * labels are omitted in 3D (Plotly 3D annotations are awkward and the Wolfram
 * default 3D graph shows none).
 *
 * Memory (SPEC section 4): returns a fresh Graphics3D tree; read-only over its
 * arguments.
 */

#include "graph.h"
#include "expr.h"
#include "sym_names.h"
#include <stdlib.h>

static const double D3_VERTEX[3] = { 0.24, 0.52, 0.90 };  /* blue             */
static const double D3_EDGE[3]   = { 0.55, 0.60, 0.72 };  /* slate            */
static const double D3_ACCENT[3] = { 0.95, 0.45, 0.10 };  /* highlight orange */
static const double D3_DIM[3]    = { 0.72, 0.75, 0.82 };  /* dimmed grey      */

static Expr* point3(double x, double y, double z) {
    Expr* xyz[3] = { expr_new_real(x), expr_new_real(y), expr_new_real(z) };
    return expr_new_function(expr_new_symbol(SYM_List), xyz, 3);
}

static Expr* rgb3(const double c[3]) {
    Expr* a[3] = { expr_new_real(c[0]), expr_new_real(c[1]), expr_new_real(c[2]) };
    return expr_new_function(expr_new_symbol(SYM_RGBColor), a, 3);
}

Expr* graph_render3d(const Expr* g, const GraphStyle* st) {
    if (!graph_is_valid_head(g, SYM_Graph3D)) return NULL;
    GraphStyle defaults = {0};
    defaults.show_labels = 0;
    if (!st) st = &defaults;

    const Expr* verts = g->data.function.args[0];
    const Expr* edges = g->data.function.args[1];
    int n = (int)verts->data.function.arg_count;
    size_t ne = edges->data.function.arg_count;

    double* x = (n > 0) ? calloc((size_t)n, sizeof(double)) : NULL;
    double* y = (n > 0) ? calloc((size_t)n, sizeof(double)) : NULL;
    double* z = (n > 0) ? calloc((size_t)n, sizeof(double)) : NULL;
    if (n > 0 && (!x || !y || !z)) { free(x); free(y); free(z); return NULL; }
    graph_compute_layout3d(g, st->layout, x, y, z);

    /* (color + line) per edge, then up to two vertex Point sets (normal +
     * highlighted), each with a color directive. */
    size_t cap = ne * 2 + 4;
    Expr** prims = (cap > 0) ? calloc(cap, sizeof(Expr*)) : NULL;
    size_t p = 0;

    for (size_t kk = 0; kk < ne; kk++) {
        const Expr* e = edges->data.function.args[kk];
        int iu = graph_vertex_index(verts, e->data.function.args[0]);
        int iv = graph_vertex_index(verts, e->data.function.args[1]);
        if (iu < 0 || iv < 0) continue;
        const double* ec = st->hi_edge ? (st->hi_edge[kk] ? D3_ACCENT : D3_DIM) : D3_EDGE;
        prims[p++] = rgb3(ec);
        Expr* pts[2] = { point3(x[iu], y[iu], z[iu]), point3(x[iv], y[iv], z[iv]) };
        Expr* seg = expr_new_function(expr_new_symbol(SYM_List), pts, 2);
        Expr* la[1] = { seg };
        prims[p++] = expr_new_function(expr_new_symbol(SYM_Line), la, 1);
    }

    /* Vertices: one Point set for the normal color, one for highlighted (if a
     * highlight mask is present). */
    for (int pass = 0; pass < 2; pass++) {
        int want_hi = (pass == 1);
        if (want_hi && !st->hi_vert) break;
        int cnt = 0;
        for (int i = 0; i < n; i++) {
            int hi = st->hi_vert ? st->hi_vert[i] : 0;
            if (st->hi_vert && (hi != want_hi)) continue;
            cnt++;
        }
        if (cnt == 0) continue;
        Expr** coords = calloc((size_t)cnt, sizeof(Expr*));
        int ci = 0;
        for (int i = 0; i < n; i++) {
            int hi = st->hi_vert ? st->hi_vert[i] : 0;
            if (st->hi_vert && (hi != want_hi)) continue;
            coords[ci++] = point3(x[i], y[i], z[i]);
        }
        const double* vc = st->hi_vert ? (want_hi ? D3_ACCENT : D3_DIM) : D3_VERTEX;
        prims[p++] = rgb3(vc);
        Expr* plist = expr_new_function(expr_new_symbol(SYM_List), coords, (size_t)cnt);
        free(coords);
        Expr* pa[1] = { plist };
        prims[p++] = expr_new_function(expr_new_symbol(SYM_Point), pa, 1);
        if (!st->hi_vert) break;   /* single uniform set */
    }

    free(x); free(y); free(z);
    Expr* prim_list = expr_new_function(expr_new_symbol(SYM_List), prims, p);
    free(prims);
    Expr* ga[1] = { prim_list };
    return expr_new_function(expr_new_symbol(SYM_Graphics3D), ga, 1);
}

Expr* graph_default_graphics3d(const Expr* g) {
    return graph_render3d(g, NULL);
}