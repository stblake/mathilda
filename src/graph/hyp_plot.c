/* hyp_plot.c - HypergraphPlot[h, opts]: draw a hypergraph as Graphics[...].
 *
 * LAYOUT. The vertices are placed by the stress layout (glayout.c) of the star
 * expansion: one extra node per hyperedge of two or more vertices, joined to
 * each of its members. Members of a hyperedge therefore sit around a common
 * centre at distance about 1 from it, and hyperedges that share vertices are
 * drawn next to each other. GraphLayout picks another embedding of the same
 * star expansion; VertexCoordinates overrides any subset of the positions.
 *
 * DRAWING. Each hyperedge is the convex hull of its members, inflated outward
 * by a margin with rounded corners (the Minkowski sum with a disk, sampled
 * every 15 degrees): a translucent filled Polygon plus an opaque outline in
 * one colour of the ColorData[97] palette, cycled by hyperedge index. So a
 * hyperedge of size 2 is a stadium (a thick translucent line with round caps)
 * and one of size 1 a circle around its vertex. Hyperedges are drawn largest
 * first so small ones stay visible on top; a hyperedge sharing a vertex with
 * earlier (smaller) ones gets a larger margin, so nested hyperedges show as
 * concentric outlines rather than coincident ones. Vertex disks and labels go
 * on top, as in GraphPlot.
 *
 * Memory (SPEC section 4): returns a freshly-allocated Graphics tree; the
 * evaluator frees res.
 */

#include "graph_hyper.h"
#include "glayout.h"
#include "expr.h"
#include "eval.h"
#include "sym_names.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define HP_FILL_OPACITY 0.22

static int is_head(const Expr* e, const char* sym) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == sym;
}

static int is_rule(const Expr* e) {
    return (is_head(e, SYM_Rule) || is_head(e, SYM_RuleDelayed))
        && e->data.function.arg_count == 2;
}

static int hyp_vindex(const void* ctx, const Expr* v) {
    return graph_vidx_get((const GraphVIdx*)ctx, v);
}

static int cmp_pt(const void* a, const void* b) {
    const double* p = (const double*)a; const double* q = (const double*)b;
    if (p[0] != q[0]) return p[0] < q[0] ? -1 : 1;
    if (p[1] != q[1]) return p[1] < q[1] ? -1 : 1;
    return 0;
}

static double cross3(const double* o, const double* a, const double* b) {
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
}

/* Convex hull (Andrew's monotone chain) of the k points in pts (sorted in
 * place), counter-clockwise, collinear points dropped, into hull (2k + 2
 * doubles). Returns the hull size (1 for coincident points). */
static int convex_hull(double* pts, int k, double* hull) {
    qsort(pts, (size_t)k, 2 * sizeof(double), cmp_pt);
    int u = 0;
    for (int i = 0; i < k; i++) {                      /* dedupe */
        if (u > 0 && fabs(pts[2 * i] - pts[2 * (u - 1)]) < 1e-12
                  && fabs(pts[2 * i + 1] - pts[2 * (u - 1) + 1]) < 1e-12) continue;
        pts[2 * u] = pts[2 * i]; pts[2 * u + 1] = pts[2 * i + 1]; u++;
    }
    k = u;
    if (k <= 2) { memcpy(hull, pts, sizeof(double) * 2 * (size_t)k); return k; }
    int h = 0;
    for (int i = 0; i < k; i++) {
        while (h >= 2 && cross3(hull + 2 * (h - 2), hull + 2 * (h - 1), pts + 2 * i) <= 1e-12) h--;
        hull[2 * h] = pts[2 * i]; hull[2 * h + 1] = pts[2 * i + 1]; h++;
    }
    for (int i = k - 2, lo = h + 1; i >= 0; i--) {
        while (h >= lo && cross3(hull + 2 * (h - 2), hull + 2 * (h - 1), pts + 2 * i) <= 1e-12) h--;
        hull[2 * h] = pts[2 * i]; hull[2 * h + 1] = pts[2 * i + 1]; h++;
    }
    return h - 1;                                       /* last == first */
}

/* The hull inflated by delta with round corners, as a closed point List;
 * extends bb. */
static Expr* rounded_hull(const double* hull, int h, double delta, double* bb, double* area) {
    enum { STEP_DEG = 15 };
    Expr** pts = malloc(sizeof(Expr*) * (size_t)(h * (360 / STEP_DEG + 2) + 4));
    if (!pts) return NULL;
    size_t np = 0;
    double a2 = 0;
    for (int i = 0; i < h; i++) {
        const double* p = hull + 2 * i;
        double a0, a1;
        if (h == 1) { a0 = 0; a1 = 2 * M_PI; }
        else {
            const double* pp = hull + 2 * ((i + h - 1) % h);
            const double* pn = hull + 2 * ((i + 1) % h);
            /* outward normals of the incoming and outgoing edges (CCW hull) */
            a0 = atan2(-(p[0] - pp[0]), p[1] - pp[1]);
            a1 = atan2(-(pn[0] - p[0]), pn[1] - p[1]);
            while (a1 < a0) a1 += 2 * M_PI;
            if (h == 2 && a1 - a0 < 1e-9) a1 += M_PI;   /* degenerate safety */
        }
        int steps = (int)ceil((a1 - a0) / (STEP_DEG * M_PI / 180.0));
        if (steps < 1) steps = 1;
        for (int s = 0; s <= steps; s++) {
            if (h == 1 && s == steps) break;
            double t = a0 + (a1 - a0) * s / steps;
            double x = p[0] + delta * cos(t), y = p[1] + delta * sin(t);
            pts[np++] = gd_pt(x, y);
            if (x < bb[0]) bb[0] = x;
            if (x > bb[1]) bb[1] = x;
            if (y < bb[2]) bb[2] = y;
            if (y > bb[3]) bb[3] = y;
        }
    }
    /* area of the inflated shape: hull area + perimeter*delta + pi delta^2 */
    double per = 0;
    for (int i = 0; i < h && h > 1; i++) {
        const double* p = hull + 2 * i; const double* q = hull + 2 * ((i + 1) % h);
        a2 += p[0] * q[1] - q[0] * p[1];
        per += sqrt((q[0] - p[0]) * (q[0] - p[0]) + (q[1] - p[1]) * (q[1] - p[1]));
    }
    *area = fabs(a2) / 2 + per * delta + M_PI * delta * delta;
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), pts, np);
    free(pts);
    return list;
}

typedef struct { int j; double area; } HOrder;

static int cmp_horder(const void* a, const void* b) {
    const HOrder* x = (const HOrder*)a; const HOrder* y = (const HOrder*)b;
    if (x->area != y->area) return x->area > y->area ? -1 : 1;
    return x->j - y->j;
}

static const char* const HP_CONSUMED[] = {
    "GraphLayout", "VertexCoordinates", "VertexLabels", "VertexStyle",
    "VertexSize", "EdgeStyle", NULL
};

Expr* builtin_hypergraph_plot(Expr* res) {
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;
    for (size_t i = 1; i < argc; i++)
        if (!is_rule(res->data.function.args[i])) return NULL;
    const Expr* h = res->data.function.args[0];
    Expr* built = NULL;
    if (is_head(h, SYM_List)) {
        Expr* a[1] = { expr_copy((Expr*)h) };
        Expr* call = expr_new_function(expr_new_symbol(hyp_sym_hypergraph()), a, 1);
        built = evaluate(call);
        expr_free(call);
        h = built;
    }
    HypView V;
    if (!hyp_view(h, &V, 0)) { if (built) expr_free(built); return NULL; }
    int n = V.n, m = V.m;

    /* ---- star expansion ------------------------------------------------ */
    int* node = malloc(sizeof(int) * (size_t)(m + 1));   /* star node per hyperedge, or -1 */
    int N = n, M = 0;
    for (int j = 0; j < m; j++) {
        int sz = V.soff[j + 1] - V.soff[j];
        node[j] = sz >= 2 ? N++ : -1;
        if (sz >= 2) M += sz;
    }
    int* eu = malloc(sizeof(int) * (size_t)(M + 1));
    int* ev = malloc(sizeof(int) * (size_t)(M + 1));
    double* xy = calloc(2 * (size_t)N + 2, sizeof(double));
    if (!node || !eu || !ev || !xy) {
        free(node); free(eu); free(ev); free(xy); if (built) expr_free(built); return NULL;
    }
    { int t = 0;
      for (int j = 0; j < m; j++) {
          if (node[j] < 0) continue;
          for (int a = V.soff[j]; a < V.soff[j + 1]; a++) { eu[t] = V.sv[a]; ev[t] = node[j]; t++; }
      } }
    GLMethod method = GL_STRESS;
    const Expr* lo = gd_option(res, 1, "GraphLayout");
    if (lo && (!glayout_parse_method(lo, &method) || method == GL_AUTOMATIC)) method = GL_STRESS;
    const Expr* vc = gd_option(res, 1, "VertexCoordinates");
    int given = vc ? gd_apply_vertex_coordinates(vc, n, xy, hyp_vindex, V.ix) : 0;
    if (given < n) {
        if (glayout_compute(N, M, eu, ev, NULL, method, xy) < 0) {
            free(node); free(eu); free(ev); free(xy); if (built) expr_free(built); return NULL;
        }
        if (vc) gd_apply_vertex_coordinates(vc, n, xy, hyp_vindex, V.ix);
    }

    /* ---- sizes: spacing from the vertex positions alone ----------------- */
    double unit = gd_drawing_unit(n, 0, NULL, NULL, xy);
    double bb[4];
    double r = gd_vertex_radius(n, xy, unit, 0, bb);
    if (n == 0) { bb[0] = bb[2] = -1; bb[1] = bb[3] = 1; }
    else { bb[0] -= 1.15 * r; bb[1] += 1.15 * r; bb[2] -= 1.15 * r; bb[3] += 1.15 * r; }
    gd_min_extent(bb, 2 * unit);
    double delta0 = 0.24 * unit, dstep = 0.10 * unit;
    if (delta0 < 2.2 * r) delta0 = 2.2 * r;

    /* Margin levels: process hyperedges by increasing size; each takes one
     * more than the most hyperedges already drawn around any of its vertices. */
    int* cnt = calloc((size_t)n + 1, sizeof(int));
    int* level = calloc((size_t)m + 1, sizeof(int));
    HOrder* ord = malloc(sizeof(HOrder) * (size_t)(m + 1));
    double* pts = malloc(sizeof(double) * 2 * (size_t)(n + 1));
    double* hull = malloc(sizeof(double) * (4 * (size_t)n + 8));   /* chain: <= 2k points */
    Expr** shapes = calloc((size_t)m + 1, sizeof(Expr*));
    GDPrims P = {0};
    int ok = cnt && level && ord && pts && hull && shapes;
    if (ok) {
        for (int j = 0; j < m; j++) { ord[j].j = j; ord[j].area = (double)(V.soff[j + 1] - V.soff[j]); }
        /* increasing size, then index */
        for (int a = 1; a < m; a++) {
            HOrder t = ord[a]; int b = a - 1;
            while (b >= 0 && (ord[b].area > t.area || (ord[b].area == t.area && ord[b].j > t.j))) {
                ord[b + 1] = ord[b]; b--;
            }
            ord[b + 1] = t;
        }
        for (int t = 0; t < m; t++) {
            int j = ord[t].j, lv = 0;
            for (int a = V.soff[j]; a < V.soff[j + 1]; a++) if (cnt[V.sv[a]] > lv) lv = cnt[V.sv[a]];
            level[j] = lv;
            for (int a = V.soff[j]; a < V.soff[j + 1]; a++) cnt[V.sv[a]]++;
        }
        /* shapes */
        for (int j = 0; j < m; j++) {
            int k = 0;
            for (int a = V.soff[j]; a < V.soff[j + 1]; a++) {
                pts[2 * k] = xy[2 * V.sv[a]]; pts[2 * k + 1] = xy[2 * V.sv[a] + 1]; k++;
            }
            ord[j].j = j; ord[j].area = 0;
            if (k == 0) continue;                          /* empty hyperedge */
            int hs = convex_hull(pts, k, hull);
            shapes[j] = rounded_hull(hull, hs, delta0 + dstep * level[j], bb, &ord[j].area);
            if (!shapes[j]) ok = 0;
        }
        qsort(ord, (size_t)m, sizeof(HOrder), cmp_horder);
    }
    if (ok) {
        double pw, ph;
        /* label avoidance: directions to the centroids of v's hyperedges */
        char** texts = gd_vertex_texts(gd_option(res, 1, "VertexLabels"), n, V.verts);
        int* aoff = NULL; double* ang = NULL;
        if (texts) {
            aoff = calloc((size_t)n + 1, sizeof(int));
            ang = malloc(sizeof(double) * (size_t)(V.soff[m] + 1));
            if (aoff && ang) {
                for (int j = 0; j < m; j++)
                    if (V.soff[j + 1] - V.soff[j] >= 2)
                        for (int a = V.soff[j]; a < V.soff[j + 1]; a++) aoff[V.sv[a] + 1]++;
                for (int i = 0; i < n; i++) aoff[i + 1] += aoff[i];
                int* fill = malloc(sizeof(int) * (size_t)(n + 1));
                if (fill) {
                    memcpy(fill, aoff, sizeof(int) * (size_t)n);
                    for (int j = 0; j < m; j++) {
                        int sz = V.soff[j + 1] - V.soff[j];
                        if (sz < 2) continue;
                        double cx = 0, cy = 0;
                        for (int a = V.soff[j]; a < V.soff[j + 1]; a++) { cx += xy[2 * V.sv[a]]; cy += xy[2 * V.sv[a] + 1]; }
                        cx /= sz; cy /= sz;
                        for (int a = V.soff[j]; a < V.soff[j + 1]; a++) {
                            int v = V.sv[a];
                            ang[fill[v]++] = atan2(cy - xy[2 * v + 1], cx - xy[2 * v]);
                        }
                    }
                    free(fill);
                } else { free(aoff); free(ang); aoff = NULL; ang = NULL; }
            } else { free(aoff); free(ang); aoff = NULL; ang = NULL; }
        }
        GDPrims labels = {0};
        gd_frame(bb, res, 1, n, xy, r, texts, aoff, ang, &labels, &pw, &ph);
        double wpt = (pw - 20.0) / (bb[1] - bb[0]), hpt = (ph - 20.0) / (bb[3] - bb[2]);
        if (hpt < wpt) wpt = hpt;
        double plot_w_pt = (bb[1] - bb[0]) * wpt;

        gd_push(&P, gd_head1(SYM_Thickness, expr_new_real(1.0 / plot_w_pt)));
        for (int t = 0; t < m; t++) {
            int j = ord[t].j;
            if (!shapes[j]) continue;
            double cr, cg, cb;
            gd_palette(j + 1, &cr, &cg, &cb);      /* skip the vertex blue */
            gd_push(&P, gd_rgb(cr, cg, cb));
            gd_push(&P, gd_head1(SYM_Opacity, expr_new_real(HP_FILL_OPACITY)));
            gd_push(&P, gd_head1(SYM_Polygon, expr_copy(shapes[j])));
            gd_push(&P, gd_head1(SYM_Opacity, expr_new_real(1.0)));
            /* outline: the same points, closed */
            Expr* ring = shapes[j];
            size_t k = ring->data.function.arg_count;
            Expr** lp = malloc(sizeof(Expr*) * (k + 1));
            if (!lp) { P.oom = 1; break; }
            for (size_t q = 0; q < k; q++) lp[q] = expr_copy(ring->data.function.args[q]);
            lp[k] = expr_copy(ring->data.function.args[0]);
            Expr* line = expr_new_function(expr_new_symbol(SYM_List), lp, k + 1);
            free(lp);
            gd_push(&P, gd_rgb(cr * 0.8, cg * 0.8, cb * 0.8));
            gd_push(&P, gd_head1(SYM_Line, line));
        }
        /* vertices */
        const Expr* vsopt = gd_option(res, 1, "VertexStyle");
        const Expr** vstyle = n ? calloc((size_t)n, sizeof(Expr*)) : NULL;
        for (int i = 0; i < n && vstyle; i++)
            vstyle[i] = gd_style_for(vsopt, V.verts->data.function.args[i], gd_vertex_match);
        Expr* dark = gd_rgb(0.25, 0.30, 0.42);
        for (int i = 0; i < n && vstyle; i++) if (!vstyle[i]) vstyle[i] = dark;
        if (n) gd_emit_vertices(&P, n, xy, r, vstyle, NULL, plot_w_pt);
        expr_free(dark);
        free((void*)vstyle);
        for (size_t i = 0; i < labels.n; i++) gd_push(&P, labels.p[i]);
        free(labels.p);
        gd_free_texts(texts, n); free(aoff); free(ang);
        Expr* out = P.oom ? NULL : gd_finish(&P, bb, pw, ph, res, 1, HP_CONSUMED);
        gd_prims_free(&P);
        for (int j = 0; j < m; j++) if (shapes[j]) expr_free(shapes[j]);
        free(shapes); free(cnt); free(level); free(ord); free(pts); free(hull);
        free(node); free(eu); free(ev); free(xy);
        if (built) expr_free(built);
        return out;
    }
    gd_prims_free(&P);
    if (shapes) for (int j = 0; j < m; j++) if (shapes[j]) expr_free(shapes[j]);
    free(shapes); free(cnt); free(level); free(ord); free(pts); free(hull);
    free(node); free(eu); free(ev); free(xy);
    if (built) expr_free(built);
    return NULL;
}
