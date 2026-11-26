/*
 * geometry_regions.c -- the region forms the polygon kernel of src/geometry.c
 * does not read directly, reduced to the ones it does.
 *
 *   Rectangle[{x1, y1}, {x2, y2}] / Rectangle[{x, y}] / Rectangle[]
 *   Triangle[{p1, p2, p3}] / Triangle[]
 *         rewritten as the equivalent Polygon (Area, Perimeter,
 *         RegionCentroid, RegionMember, RegionMeasure).
 *   Polygon[{poly1, poly2, ...}]          a multi-polygon: measures add, the
 *                                         centroid is area-weighted, a point
 *                                         is a member of any part.
 *   Polygon[outer -> {hole1, ...}]        a polygon with holes: the holes'
 *                                         areas subtract, their boundaries
 *                                         add to the perimeter, a point is a
 *                                         member when it is in the outer
 *                                         polygon and not strictly inside a
 *                                         hole (the hole's edge belongs to
 *                                         the region).
 *   Polygon[{{x, y, z}, ...}]             RegionMember of a planar 3-D
 *                                         polygon: coplanar, and inside the
 *                                         projection onto the coordinate
 *                                         plane the normal is most along.
 *   RegionMember on Point, Line (any dimension), Circle[c, r] and Ball[c, r]
 *   (any dimension).
 *
 * Exact coordinates that are not rational (Sqrt[2], Pi) take the symbolic
 * path below: the same shoelace / crossing / monotone-chain algorithms, with
 * the arithmetic built as expressions and every branch decided by the
 * certified Sign[] of an exact numeric quantity.  Inexact coordinates mixed
 * with such constants are numericized first (the machine path then runs).
 *
 * Structural region code, not an element-wise numeric map: no NDArray or
 * Compile surface (EXEMPT, as for the rest of the geometry heads).
 */

#include "geometry.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"       /* head_is */
#include "eval.h"
#include "expr.h"
#include "ndarray.h"
#include "numeric.h"      /* expr_is_numeric_quantity */
#include "sym_names.h"

/* ---- expression helpers (each consumes its arguments) ------------------ */

static Expr* gr_fn(const char* h, Expr** a, size_t n) {
    return expr_new_function(expr_new_symbol(h), a, n);
}
static Expr* gr_fn1(const char* h, Expr* a) { return gr_fn(h, (Expr*[]){ a }, 1); }
static Expr* gr_fn2(const char* h, Expr* a, Expr* b) { return gr_fn(h, (Expr*[]){ a, b }, 2); }
static Expr* gr_c(const Expr* e) { return expr_copy((Expr*)e); }
static Expr* gr_ev(Expr* e) { return eval_and_free(e); }
static Expr* gr_add(Expr* a, Expr* b) { return gr_fn2(SYM_Plus, a, b); }
static Expr* gr_mul(Expr* a, Expr* b) { return gr_fn2(SYM_Times, a, b); }
static Expr* gr_sub(Expr* a, Expr* b) { return gr_add(a, gr_mul(expr_new_integer(-1), b)); }
static Expr* gr_inv(Expr* a) { return gr_fn2(SYM_Power, a, expr_new_integer(-1)); }
static size_t gr_len(const Expr* e) { return e->data.function.arg_count; }
static const Expr* gr_at(const Expr* e, size_t i) { return e->data.function.args[i]; }

static Expr* gr_list(Expr** a, size_t n) { return gr_fn(SYM_List, a, n); }

/* 1 True, 0 False, -1 undecided; consumes e. */
static int gr_truth(Expr* e) {
    int t = -1;
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol.name == SYM_True) t = 1;
        else if (e->data.symbol.name == SYM_False) t = 0;
    }
    expr_free(e);
    return t;
}

static bool gr_has_inexact(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < gr_len(e); i++) if (gr_has_inexact(gr_at(e, i))) return true;
    return false;
}

static bool gr_is_rational_leaf(const Expr* e) {
    if (expr_is_integer_like(e)) return true;
    return head_is(e, SYM_Rational);
}

/* A point: a List of `d` (0: any) numeric quantities. */
static bool gr_is_point(const Expr* e, size_t d) {
    if (!head_is(e, SYM_List) || gr_len(e) == 0) return false;
    if (d && gr_len(e) != d) return false;
    for (size_t i = 0; i < gr_len(e); i++) {
        const Expr* x = gr_at(e, i);
        if (head_is(x, SYM_List) || !expr_is_numeric_quantity(x)) return false;
    }
    return true;
}

/* A nonempty List of points of one dimension; *d receives it. */
static bool gr_is_point_list(const Expr* e, size_t* d) {
    if (!head_is(e, SYM_List) || gr_len(e) == 0) return false;
    if (!head_is(gr_at(e, 0), SYM_List)) return false;
    size_t k = gr_len(gr_at(e, 0));
    for (size_t i = 0; i < gr_len(e); i++) if (!gr_is_point(gr_at(e, i), k)) return false;
    if (d) *d = k;
    return true;
}

/* Every element a point list (a multi-polygon / several lines). */
static bool gr_is_multi(const Expr* e) {
    if (!head_is(e, SYM_List) || gr_len(e) == 0) return false;
    for (size_t i = 0; i < gr_len(e); i++) {
        size_t d;
        if (!gr_is_point_list(gr_at(e, i), &d)) return false;
    }
    return true;
}

static Expr* gr_unpack(const Expr* e) {
    if (is_ndarray(e)) {
        Expr* m = ndarray_to_nested_list(e);
        if (m) return m;
    }
    return gr_c(e);
}

/* h[region, extra] evaluated; NULL when it stays unevaluated (head h). */
static Expr* gr_call(const char* h, Expr* region, const Expr* extra) {
    Expr* call = extra ? gr_fn2(h, region, gr_c(extra)) : gr_fn1(h, region);
    Expr* r = gr_ev(call);
    if (head_is(r, h)) { expr_free(r); return NULL; }
    return r;
}

/* ---- Rectangle / Triangle as Polygon ----------------------------------- */

static Expr* gr_pt2(Expr* x, Expr* y) { return gr_list((Expr*[]){ x, y }, 2); }

static Expr* gr_as_polygon(const Expr* g) {
    if (head_is(g, SYM_Rectangle)) {
        size_t n = gr_len(g);
        Expr *x1, *y1, *x2, *y2;
        if (n == 0) {
            x1 = expr_new_integer(0); y1 = expr_new_integer(0);
            x2 = expr_new_integer(1); y2 = expr_new_integer(1);
        } else if (n <= 2) {
            const Expr* a = gr_at(g, 0);
            if (!head_is(a, SYM_List) || gr_len(a) != 2) return NULL;
            x1 = gr_c(gr_at(a, 0)); y1 = gr_c(gr_at(a, 1));
            if (n == 2) {
                const Expr* b = gr_at(g, 1);
                if (!head_is(b, SYM_List) || gr_len(b) != 2) { expr_free(x1); expr_free(y1); return NULL; }
                x2 = gr_c(gr_at(b, 0)); y2 = gr_c(gr_at(b, 1));
            } else {
                x2 = gr_ev(gr_add(gr_c(x1), expr_new_integer(1)));
                y2 = gr_ev(gr_add(gr_c(y1), expr_new_integer(1)));
            }
        } else {
            return NULL;
        }
        Expr* pts = gr_list((Expr*[]){ gr_pt2(gr_c(x1), gr_c(y1)), gr_pt2(gr_c(x2), gr_c(y1)),
                                       gr_pt2(gr_c(x2), gr_c(y2)), gr_pt2(gr_c(x1), gr_c(y2)) }, 4);
        expr_free(x1); expr_free(y1); expr_free(x2); expr_free(y2);
        return gr_fn1(SYM_Polygon, pts);
    }
    if (head_is(g, SYM_Triangle)) {
        if (gr_len(g) == 0) {
            Expr* pts = gr_list((Expr*[]){ gr_pt2(expr_new_integer(0), expr_new_integer(0)),
                                           gr_pt2(expr_new_integer(1), expr_new_integer(0)),
                                           gr_pt2(expr_new_integer(0), expr_new_integer(1)) }, 3);
            return gr_fn1(SYM_Polygon, pts);
        }
        size_t d;
        if (gr_len(g) == 1 && gr_is_point_list(gr_at(g, 0), &d) && gr_len(gr_at(g, 0)) == 3)
            return gr_fn1(SYM_Polygon, gr_c(gr_at(g, 0)));
    }
    return NULL;
}

/* ---- symbolic exact polygons (Sqrt[2], Pi coordinates) ------------------ */

/* Certified sign of an exact numeric quantity: -1, 0, 1, or 2 when Sign does
 * not decide it.  Consumes e. */
static int gr_sign(Expr* e) {
    Expr* s = gr_ev(gr_fn1(SYM_Sign, gr_ev(gr_fn1(SYM_Expand, e))));
    int r = 2;
    if (s->type == EXPR_INTEGER && s->data.integer >= -1 && s->data.integer <= 1)
        r = (int)s->data.integer;
    expr_free(s);
    return r;
}

/* Polygon vertex list fit for the symbolic path: 2-D, every leaf an exact
 * numeric quantity. */
static bool gr_exact_points(const Expr* pts) {
    size_t d;
    return gr_is_point_list(pts, &d) && d == 2 && !gr_has_inexact(pts);
}

static bool gr_all_rational(const Expr* e) {
    if (e->type == EXPR_FUNCTION && head_is(e, SYM_List)) {
        for (size_t i = 0; i < gr_len(e); i++) if (!gr_all_rational(gr_at(e, i))) return false;
        return true;
    }
    return gr_is_rational_leaf(e);
}

#define GX(p, i) gr_at(gr_at((p), (i)), 0)
#define GY(p, i) gr_at(gr_at((p), (i)), 1)

/* Twice the signed area, expanded. */
static Expr* gr_area2(const Expr* p) {
    size_t n = gr_len(p);
    Expr** t = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        t[i] = gr_sub(gr_mul(gr_c(GX(p, i)), gr_c(GY(p, j))), gr_mul(gr_c(GX(p, j)), gr_c(GY(p, i))));
    }
    Expr* s = gr_ev(gr_fn1(SYM_Expand, gr_fn(SYM_Plus, t, n)));
    free(t);
    return s;
}

static Expr* gr_sym_area(const Expr* p) {
    if (gr_len(p) < 3) return expr_new_symbol(SYM_Undefined);
    Expr* s = gr_area2(p);
    int sg = gr_sign(gr_c(s));
    if (sg == 2) { expr_free(s); return NULL; }
    Expr* half = gr_fn2(SYM_Rational, expr_new_integer(sg), expr_new_integer(2));
    return gr_ev(gr_mul(half, s));
}

static Expr* gr_sym_perimeter(const Expr* p) {
    size_t n = gr_len(p);
    if (n < 3) return expr_new_symbol(SYM_Undefined);
    Expr** t = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        Expr* dx = gr_sub(gr_c(GX(p, j)), gr_c(GX(p, i)));
        Expr* dy = gr_sub(gr_c(GY(p, j)), gr_c(GY(p, i)));
        Expr* q = gr_add(gr_fn2(SYM_Power, dx, expr_new_integer(2)),
                         gr_fn2(SYM_Power, dy, expr_new_integer(2)));
        t[i] = gr_fn1(SYM_Sqrt, gr_fn1(SYM_Expand, q));
    }
    Expr* s = gr_ev(gr_fn(SYM_Plus, t, n));
    free(t);
    return s;
}

static Expr* gr_sym_centroid(const Expr* p) {
    size_t n = gr_len(p);
    if (n < 3) return NULL;
    Expr* a2 = gr_area2(p);
    if (gr_sign(gr_c(a2)) == 0) { expr_free(a2); return NULL; }
    Expr** sx = malloc(sizeof(Expr*) * n);
    Expr** sy = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        Expr* cr = gr_sub(gr_mul(gr_c(GX(p, i)), gr_c(GY(p, j))), gr_mul(gr_c(GX(p, j)), gr_c(GY(p, i))));
        sx[i] = gr_mul(gr_add(gr_c(GX(p, i)), gr_c(GX(p, j))), gr_c(cr));
        sy[i] = gr_mul(gr_add(gr_c(GY(p, i)), gr_c(GY(p, j))), cr);
    }
    Expr* den = gr_mul(expr_new_integer(3), a2);
    Expr* cx = gr_ev(gr_fn1(SYM_Simplify, gr_mul(gr_fn(SYM_Plus, sx, n), gr_inv(gr_c(den)))));
    Expr* cy = gr_ev(gr_fn1(SYM_Simplify, gr_mul(gr_fn(SYM_Plus, sy, n), gr_inv(den))));
    free(sx); free(sy);
    return gr_list((Expr*[]){ cx, cy }, 2);
}

/* sign of (b - a) x (q - a) */
static int gr_cross_sign(const Expr* ax, const Expr* ay, const Expr* bx, const Expr* by,
                         const Expr* qx, const Expr* qy) {
    Expr* u = gr_mul(gr_sub(gr_c(bx), gr_c(ax)), gr_sub(gr_c(qy), gr_c(ay)));
    Expr* v = gr_mul(gr_sub(gr_c(by), gr_c(ay)), gr_sub(gr_c(qx), gr_c(ax)));
    return gr_sign(gr_sub(u, v));
}

static int gr_cmp(const Expr* a, const Expr* b) { return gr_sign(gr_sub(gr_c(a), gr_c(b))); }

/* Boundary-inclusive even-odd membership; -1 undecided. */
static int gr_sym_member(const Expr* p, const Expr* q) {
    size_t n = gr_len(p);
    const Expr* qx = gr_at(q, 0);
    const Expr* qy = gr_at(q, 1);
    int inside = 0;
    for (size_t i = 0; i < n; i++) {
        size_t j = (i + 1) % n;
        int cs = gr_cross_sign(GX(p, i), GY(p, i), GX(p, j), GY(p, j), qx, qy);
        if (cs == 2) return -1;
        if (cs == 0) {
            int a = gr_cmp(GX(p, i), qx), b = gr_cmp(qx, GX(p, j));
            int c = gr_cmp(GX(p, j), qx), d = gr_cmp(qx, GX(p, i));
            int e = gr_cmp(GY(p, i), qy), f = gr_cmp(qy, GY(p, j));
            int g = gr_cmp(GY(p, j), qy), h = gr_cmp(qy, GY(p, i));
            if (a == 2 || b == 2 || c == 2 || d == 2 || e == 2 || f == 2 || g == 2 || h == 2) return -1;
            bool in_x = (a <= 0 && b <= 0) || (c <= 0 && d <= 0);
            bool in_y = (e <= 0 && f <= 0) || (g <= 0 && h <= 0);
            if (in_x && in_y) return 1;
        }
        int ai = gr_cmp(GY(p, i), qy), bj = gr_cmp(GY(p, j), qy);
        if (ai == 2 || bj == 2) return -1;
        if ((ai > 0) != (bj > 0)) {
            int s = gr_cmp(GY(p, j), GY(p, i));
            if (s == 2) return -1;
            if (cs == s) inside = !inside;
        }
    }
    return inside;
}

/* Monotone chain over exact numeric points (ConvexHullRegion). */
static const Expr* g_hull_pts;
static int g_hull_undecided;

static int gr_idx_cmp(const void* a, const void* b) {
    size_t i = *(const size_t*)a, j = *(const size_t*)b;
    int c = gr_cmp(GX(g_hull_pts, i), GX(g_hull_pts, j));
    if (c == 2) { g_hull_undecided = 1; return 0; }
    if (c) return c;
    c = gr_cmp(GY(g_hull_pts, i), GY(g_hull_pts, j));
    if (c == 2) { g_hull_undecided = 1; return 0; }
    return c;
}

static int gr_turn(const Expr* p, size_t a, size_t b, size_t c) {
    int s = gr_cross_sign(GX(p, a), GY(p, a), GX(p, b), GY(p, b), GX(p, c), GY(p, c));
    if (s == 2) g_hull_undecided = 1;
    return s == 2 ? 0 : s;
}

Expr* geom_symbolic_hull(const Expr* pts) {
    if (!gr_exact_points(pts) || gr_all_rational(pts)) return NULL;
    size_t n = gr_len(pts);
    size_t* idx = malloc(n * sizeof(size_t));
    size_t* hull = malloc((2 * n + 1) * sizeof(size_t));
    for (size_t i = 0; i < n; i++) idx[i] = i;
    g_hull_pts = pts;
    g_hull_undecided = 0;
    qsort(idx, n, sizeof(size_t), gr_idx_cmp);
    size_t m = 0;
    for (size_t i = 0; i < n; i++)
        if (m == 0 || gr_idx_cmp(&idx[m - 1], &idx[i]) != 0) idx[m++] = idx[i];
    size_t k = 0;
    if (m == 1) {
        hull[k++] = idx[0];
    } else {
        for (size_t i = 0; i < m; i++) {
            while (k >= 2 && gr_turn(pts, hull[k - 2], hull[k - 1], idx[i]) <= 0) k--;
            hull[k++] = idx[i];
        }
        size_t lower = k + 1;
        for (size_t ii = m - 1; ii-- > 0; ) {
            while (k >= lower && gr_turn(pts, hull[k - 2], hull[k - 1], idx[ii]) <= 0) k--;
            hull[k++] = idx[ii];
        }
        k--;
    }
    g_hull_pts = NULL;
    Expr* out = NULL;
    if (!g_hull_undecided) {
        Expr** v = malloc(sizeof(Expr*) * (k ? k : 1));
        for (size_t i = 0; i < k; i++) v[i] = gr_c(gr_at(pts, hull[i]));
        if (k == 1) out = gr_fn(SYM_Point, v, 1);
        else out = gr_fn1(k == 2 ? SYM_Line : SYM_Polygon, gr_list(v, k));
        free(v);
    }
    free(idx); free(hull);
    return out;
}

/* ---- RegionMember on the further regions ------------------------------- */

/* Squared distance |a - b|^2 (same-length point Lists). */
static Expr* gr_dist2(const Expr* a, const Expr* b) {
    size_t d = gr_len(a);
    Expr** t = malloc(sizeof(Expr*) * d);
    for (size_t i = 0; i < d; i++)
        t[i] = gr_fn2(SYM_Power, gr_sub(gr_c(gr_at(a, i)), gr_c(gr_at(b, i))), expr_new_integer(2));
    Expr* s = gr_fn(SYM_Plus, t, d);
    free(t);
    return s;
}

/* q on the segment a-b: q - a = t (b - a), 0 <= t <= 1. -1 undecided. */
static int gr_on_segment(const Expr* a, const Expr* b, const Expr* q) {
    size_t d = gr_len(a);
    Expr* ab2 = gr_ev(gr_dist2(b, a));
    Expr** dot = malloc(sizeof(Expr*) * d);
    for (size_t i = 0; i < d; i++)
        dot[i] = gr_mul(gr_sub(gr_c(gr_at(q, i)), gr_c(gr_at(a, i))),
                        gr_sub(gr_c(gr_at(b, i)), gr_c(gr_at(a, i))));
    Expr* dt = gr_ev(gr_fn(SYM_Plus, dot, d));
    free(dot);
    int z = gr_truth(gr_ev(gr_fn2(SYM_Equal, gr_c(ab2), expr_new_integer(0))));
    if (z == 1) {           /* degenerate segment: the point itself */
        expr_free(ab2); expr_free(dt);
        return gr_truth(gr_ev(gr_fn2(SYM_Equal, gr_c(q), gr_c(a))));
    }
    if (z != 0) { expr_free(ab2); expr_free(dt); return -1; }
    Expr* t = gr_ev(gr_mul(dt, gr_inv(gr_c(ab2))));
    expr_free(ab2);
    int in = gr_truth(gr_ev(gr_fn(SYM_LessEqual, (Expr*[]){ expr_new_integer(0), gr_c(t), expr_new_integer(1) }, 3)));
    if (in != 1) { expr_free(t); return in; }
    /* the foot of the perpendicular is q itself */
    Expr** foot = malloc(sizeof(Expr*) * d);
    for (size_t i = 0; i < d; i++)
        foot[i] = gr_add(gr_c(gr_at(a, i)), gr_mul(gr_c(t), gr_sub(gr_c(gr_at(b, i)), gr_c(gr_at(a, i)))));
    Expr* f = gr_ev(gr_list(foot, d));
    free(foot);
    expr_free(t);
    int r = gr_truth(gr_ev(gr_fn2(SYM_Equal, f, gr_c(q))));
    return r;
}

/* 1 / 0 / -1 (undecided or not handled). */
static int gr_member_point(const Expr* g, const Expr* q);

static int gr_member_line(const Expr* pts, const Expr* q) {
    size_t n = gr_len(pts);
    if (n == 1) return gr_truth(gr_ev(gr_fn2(SYM_Equal, gr_c(gr_at(pts, 0)), gr_c(q))));
    int undecided = 0;
    for (size_t i = 0; i + 1 < n; i++) {
        int r = gr_on_segment(gr_at(pts, i), gr_at(pts, i + 1), q);
        if (r == 1) return 1;
        if (r < 0) undecided = 1;
    }
    return undecided ? -1 : 0;
}

/* Planar 3-D polygon: coplanar with the first non-degenerate vertex triple,
 * then the 2-D test on the projection dropping the normal's largest axis. */
static int gr_member_poly3(const Expr* pts, const Expr* q) {
    size_t n = gr_len(pts);
    if (n < 3) return -1;
    const Expr* p0 = gr_at(pts, 0);
    for (size_t k = 2; k < n; k++) {
        const Expr* p1 = gr_at(pts, 1);
        const Expr* p2 = gr_at(pts, k);
        Expr* u[3]; Expr* v[3];
        for (int i = 0; i < 3; i++) {
            u[i] = gr_sub(gr_c(gr_at(p1, i)), gr_c(gr_at(p0, i)));
            v[i] = gr_sub(gr_c(gr_at(p2, i)), gr_c(gr_at(p0, i)));
        }
        Expr* nrm[3];
        for (int i = 0; i < 3; i++) {
            int a = (i + 1) % 3, b = (i + 2) % 3;
            nrm[i] = gr_ev(gr_sub(gr_mul(gr_c(u[a]), gr_c(v[b])), gr_mul(gr_c(u[b]), gr_c(v[a]))));
        }
        for (int i = 0; i < 3; i++) { expr_free(u[i]); expr_free(v[i]); }
        Expr* nv = gr_list((Expr*[]){ gr_c(nrm[0]), gr_c(nrm[1]), gr_c(nrm[2]) }, 3);
        int zero = gr_truth(gr_ev(gr_fn2(SYM_Equal, nv, gr_list((Expr*[]){ expr_new_integer(0),
                                expr_new_integer(0), expr_new_integer(0) }, 3))));
        if (zero != 0) {
            for (int i = 0; i < 3; i++) expr_free(nrm[i]);
            if (zero < 0) return -1;
            continue;                               /* collinear triple */
        }
        /* coplanar: (q - p0) . n == 0 */
        Expr* dt = gr_add(gr_add(gr_mul(gr_sub(gr_c(gr_at(q, 0)), gr_c(gr_at(p0, 0))), gr_c(nrm[0])),
                                 gr_mul(gr_sub(gr_c(gr_at(q, 1)), gr_c(gr_at(p0, 1))), gr_c(nrm[1]))),
                          gr_mul(gr_sub(gr_c(gr_at(q, 2)), gr_c(gr_at(p0, 2))), gr_c(nrm[2])));
        int cop = gr_truth(gr_ev(gr_fn2(SYM_Equal, dt, expr_new_integer(0))));
        /* the axis the normal is most along */
        int drop = 0;
        for (int i = 1; i < 3; i++) {
            int t = gr_truth(gr_ev(gr_fn2(SYM_Greater, gr_fn1(SYM_Abs, gr_c(nrm[i])),
                                                     gr_fn1(SYM_Abs, gr_c(nrm[drop])))));
            if (t == 1) drop = i;
        }
        for (int i = 0; i < 3; i++) expr_free(nrm[i]);
        if (cop != 1) return cop;
        int a = drop == 0 ? 1 : 0, b = drop == 2 ? 1 : 2;
        Expr** pp = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++)
            pp[i] = gr_pt2(gr_c(gr_at(gr_at(pts, i), a)), gr_c(gr_at(gr_at(pts, i), b)));
        Expr* poly = gr_fn1(SYM_Polygon, gr_list(pp, n));
        free(pp);
        Expr* q2 = gr_pt2(gr_c(gr_at(q, a)), gr_c(gr_at(q, b)));
        Expr* r = gr_call(SYM_RegionMember, poly, q2);
        expr_free(q2);
        return r ? gr_truth(r) : -1;
    }
    return -1;
}

static int gr_member_point(const Expr* g, const Expr* q) {
    size_t d = gr_len(q);
    if (head_is(g, SYM_Point) && gr_len(g) == 1) {
        const Expr* a = gr_at(g, 0);
        size_t k;
        if (gr_is_point(a, d)) return gr_truth(gr_ev(gr_fn2(SYM_Equal, gr_c(a), gr_c(q))));
        if (gr_is_point_list(a, &k) && k == d) {
            int undecided = 0;
            for (size_t i = 0; i < gr_len(a); i++) {
                int r = gr_truth(gr_ev(gr_fn2(SYM_Equal, gr_c(gr_at(a, i)), gr_c(q))));
                if (r == 1) return 1;
                if (r < 0) undecided = 1;
            }
            return undecided ? -1 : 0;
        }
        return -1;
    }
    if (head_is(g, SYM_Line) && gr_len(g) == 1) {
        const Expr* a = gr_at(g, 0);
        size_t k;
        if (gr_is_point_list(a, &k) && k == d) return gr_member_line(a, q);
        if (gr_is_multi(a)) {
            int undecided = 0;
            for (size_t i = 0; i < gr_len(a); i++) {
                if (gr_len(gr_at(gr_at(a, i), 0)) != d) return -1;
                int r = gr_member_line(gr_at(a, i), q);
                if (r == 1) return 1;
                if (r < 0) undecided = 1;
            }
            return undecided ? -1 : 0;
        }
        return -1;
    }
    /* Circle[c, r] (a full circle; arcs and ellipses are not read here) and
     * Ball[c, r]: |q - c|^2 == r^2 and <= r^2. */
    bool circ = head_is(g, SYM_Circle), ball = head_is(g, SYM_Ball);
    if ((circ || ball) && gr_len(g) <= 2) {
        Expr* c;
        if (gr_len(g) >= 1) {
            if (!gr_is_point(gr_at(g, 0), d)) return -1;
            c = gr_c(gr_at(g, 0));
        } else {
            if (circ && d != 2) return -1;
            if (ball && d != 3) return -1;
            Expr** z = malloc(sizeof(Expr*) * d);
            for (size_t i = 0; i < d; i++) z[i] = expr_new_integer(0);
            c = gr_list(z, d);
            free(z);
        }
        if (circ && d != 2) { expr_free(c); return -1; }
        Expr* r = gr_len(g) == 2 ? gr_c(gr_at(g, 1)) : expr_new_integer(1);
        if (head_is(r, SYM_List) || !expr_is_numeric_quantity(r)) { expr_free(c); expr_free(r); return -1; }
        Expr* d2 = gr_dist2(q, c);
        expr_free(c);
        Expr* r2 = gr_fn2(SYM_Power, r, expr_new_integer(2));
        return gr_truth(gr_ev(gr_fn2(circ ? SYM_Equal : SYM_LessEqual, d2, r2)));
    }
    if (head_is(g, SYM_Polygon) && gr_len(g) == 1) {
        size_t k;
        if (gr_is_point_list(gr_at(g, 0), &k) && k == 3 && d == 3) return gr_member_poly3(gr_at(g, 0), q);
    }
    return -1;
}

/* ---- RegionCentroid of points, curves and solids ----------------------- */

/* sum(w_i p_i) / sum(w_i) over k points of dimension d (p[i] a List of d
 * coordinates), each coordinate Simplify'd; NULL when the weights sum to
 * zero.  w (NULL: equal weights) is borrowed. */
static Expr* gr_weighted_mean(const Expr* const* p, Expr* const* w, size_t k, size_t d) {
    Expr* wsum;
    if (w) {
        Expr** ws = malloc(sizeof(Expr*) * k);
        for (size_t i = 0; i < k; i++) ws[i] = gr_c(w[i]);
        wsum = gr_ev(gr_fn(SYM_Plus, ws, k));
        free(ws);
    } else {
        wsum = expr_new_integer((int64_t)k);
    }
    int z = gr_truth(gr_ev(gr_fn2(SYM_Equal, gr_c(wsum), expr_new_integer(0))));
    if (z != 0) { expr_free(wsum); return NULL; }
    Expr** coord = malloc(sizeof(Expr*) * d);
    Expr** t = malloc(sizeof(Expr*) * k);
    for (size_t j = 0; j < d; j++) {
        for (size_t i = 0; i < k; i++) {
            Expr* x = gr_c(gr_at(p[i], j));
            t[i] = w ? gr_mul(gr_c(w[i]), x) : x;
        }
        Expr* num = gr_fn(SYM_Plus, t, k);
        coord[j] = gr_ev(gr_fn1(SYM_Simplify, gr_mul(num, gr_inv(gr_c(wsum)))));
    }
    free(t);
    expr_free(wsum);
    Expr* out = gr_list(coord, d);
    free(coord);
    return out;
}

/* (a + b) / 2 coordinate-wise (same-length point Lists), unevaluated. */
static Expr* gr_midpoint(const Expr* a, const Expr* b) {
    size_t d = gr_len(a);
    Expr** m = malloc(sizeof(Expr*) * d);
    for (size_t i = 0; i < d; i++)
        m[i] = gr_mul(gr_fn2(SYM_Rational, expr_new_integer(1), expr_new_integer(2)),
                      gr_add(gr_c(gr_at(a, i)), gr_c(gr_at(b, i))));
    Expr* out = gr_ev(gr_list(m, d));
    free(m);
    return out;
}

/* Point[p] -> p; Point[{p1, ..., pn}] -> their mean (audit Point#2). */
static Expr* gr_centroid_point(const Expr* g) {
    if (gr_len(g) != 1) return NULL;
    const Expr* p = gr_at(g, 0);
    if (gr_is_point(p, 0)) return gr_c(p);
    size_t d;
    if (!gr_is_point_list(p, &d)) return NULL;
    return gr_weighted_mean((const Expr* const*)p->data.function.args, NULL, gr_len(p), d);
}

/* Line[{p1, ..., pn}]: the length-weighted mean of the segment midpoints
 * (audit Line#2 / RegionCentroid#4); a Line of one point is that point. */
static Expr* gr_centroid_line(const Expr* g) {
    if (gr_len(g) != 1) return NULL;
    const Expr* p = gr_at(g, 0);
    size_t d;
    if (!gr_is_point_list(p, &d)) return NULL;
    size_t n = gr_len(p);
    if (n == 1) return gr_c(gr_at(p, 0));
    Expr** mid = malloc(sizeof(Expr*) * (n - 1));
    Expr** len = malloc(sizeof(Expr*) * (n - 1));
    for (size_t i = 0; i + 1 < n; i++) {
        mid[i] = gr_midpoint(gr_at(p, i), gr_at(p, i + 1));
        len[i] = gr_ev(gr_fn1(SYM_Sqrt, gr_dist2(gr_at(p, i + 1), gr_at(p, i))));
    }
    Expr* out = gr_weighted_mean((const Expr* const*)mid, len, n - 1, d);
    for (size_t i = 0; i + 1 < n; i++) { expr_free(mid[i]); expr_free(len[i]); }
    free(mid); free(len);
    return out;
}

/* Circle[c, r] (and the ellipse Circle[c, {a, b}]): c; the circular arc
 * Circle[c, r, {t1, t2}]: c + r Sin[h]/h {Cos[m], Sin[m]}, h the half-span
 * and m the mid-angle (audit Circle#3).  The elliptical arc is left alone. */
static Expr* gr_centroid_circle(const Expr* g) {
    size_t n = gr_len(g);
    if (n == 0) return gr_pt2(expr_new_integer(0), expr_new_integer(0));
    const Expr* c = gr_at(g, 0);
    if (!gr_is_point(c, 2) && !(head_is(c, SYM_List) && gr_len(c) == 2)) return NULL;
    if (n <= 2) return gr_c(c);
    if (n != 3) return NULL;
    const Expr* r = gr_at(g, 1);
    const Expr* th = gr_at(g, 2);
    if (head_is(r, SYM_List) || !head_is(th, SYM_List) || gr_len(th) != 2) return NULL;
    Expr* span = gr_ev(gr_fn1(SYM_Abs, gr_sub(gr_c(gr_at(th, 1)), gr_c(gr_at(th, 0)))));
    if (expr_is_numeric_quantity(span)) {
        Expr* two_pi = gr_mul(expr_new_integer(2), expr_new_symbol(SYM_Pi));
        int full = gr_truth(gr_ev(gr_fn2(SYM_GreaterEqual, gr_c(span), two_pi)));
        if (full == 1) { expr_free(span); return gr_c(c); }
    }
    Expr* half = gr_fn2(SYM_Rational, expr_new_integer(1), expr_new_integer(2));
    Expr* h = gr_ev(gr_mul(gr_c(half), span));
    int zero = gr_truth(gr_ev(gr_fn2(SYM_Equal, gr_c(h), expr_new_integer(0))));
    if (zero != 0) { expr_free(h); expr_free(half); return NULL; }
    Expr* m = gr_ev(gr_mul(half, gr_add(gr_c(gr_at(th, 0)), gr_c(gr_at(th, 1)))));
    Expr* k = gr_mul(gr_ev(gr_fn1(SYM_Abs, gr_c(r))), gr_mul(gr_fn1(SYM_Sin, gr_c(h)), gr_inv(h)));
    Expr* ke = gr_ev(k);
    Expr* x = gr_add(gr_c(gr_at(c, 0)), gr_mul(gr_c(ke), gr_fn1(SYM_Cos, gr_c(m))));
    Expr* y = gr_add(gr_c(gr_at(c, 1)), gr_mul(ke, gr_fn1(SYM_Sin, m)));
    return gr_pt2(gr_ev(gr_fn1(SYM_Simplify, x)), gr_ev(gr_fn1(SYM_Simplify, y)));
}

/* Ball[c, r] / Sphere[c, r]: c (Ball[] the origin of 3-space); Cuboid[pmin]
 * the unit cube's centre pmin + 1/2, Cuboid[pmin, pmax] the midpoint
 * (audit RegionCentroid#5). */
static Expr* gr_centroid_solid(const Expr* g) {
    size_t n = gr_len(g);
    bool cub = head_is(g, SYM_Cuboid);
    if (n == 0) {
        Expr* v = cub ? gr_fn2(SYM_Rational, expr_new_integer(1), expr_new_integer(2))
                      : expr_new_integer(0);
        Expr* out = gr_list((Expr*[]){ gr_c(v), gr_c(v), gr_c(v) }, 3);
        expr_free(v);
        return out;
    }
    const Expr* c = gr_at(g, 0);
    if (!gr_is_point(c, 0) || gr_len(c) < 2) return NULL;
    if (!cub) return n <= 2 ? gr_c(c) : NULL;
    if (n == 1) {
        size_t d = gr_len(c);
        Expr** m = malloc(sizeof(Expr*) * d);
        for (size_t i = 0; i < d; i++)
            m[i] = gr_add(gr_c(gr_at(c, i)), gr_fn2(SYM_Rational, expr_new_integer(1), expr_new_integer(2)));
        Expr* out = gr_ev(gr_list(m, d));
        free(m);
        return out;
    }
    if (n != 2 || !gr_is_point(gr_at(g, 1), gr_len(c))) return NULL;
    return gr_midpoint(c, gr_at(g, 1));
}

/* A planar polygon in 3-space: the fan p0, p_i, p_{i+1} triangles weighted
 * by their area along the polygon's normal N = sum of the fan's cross
 * products (signed, so a concave polygon's reflex triangles subtract), which
 * needs no square root (audit RegionCentroid#5). */
static Expr* gr_centroid_poly3(const Expr* g) {
    if (!head_is(g, SYM_Polygon) || gr_len(g) != 1) return NULL;
    const Expr* p = gr_at(g, 0);
    size_t d;
    if (!gr_is_point_list(p, &d) || d != 3 || gr_len(p) < 3) return NULL;
    size_t n = gr_len(p);
    size_t nt = n - 2;
    Expr** cr = malloc(sizeof(Expr*) * nt * 3);     /* cross products */
    const Expr* p0 = gr_at(p, 0);
    for (size_t i = 0; i < nt; i++) {
        const Expr* a = gr_at(p, i + 1);
        const Expr* b = gr_at(p, i + 2);
        Expr* u[3]; Expr* v[3];
        for (size_t j = 0; j < 3; j++) {
            u[j] = gr_sub(gr_c(gr_at(a, j)), gr_c(gr_at(p0, j)));
            v[j] = gr_sub(gr_c(gr_at(b, j)), gr_c(gr_at(p0, j)));
        }
        for (size_t j = 0; j < 3; j++) {
            size_t j1 = (j + 1) % 3, j2 = (j + 2) % 3;
            cr[3 * i + j] = gr_ev(gr_sub(gr_mul(gr_c(u[j1]), gr_c(v[j2])), gr_mul(gr_c(u[j2]), gr_c(v[j1]))));
        }
        for (size_t j = 0; j < 3; j++) { expr_free(u[j]); expr_free(v[j]); }
    }
    Expr* N[3];
    for (size_t j = 0; j < 3; j++) {
        Expr** t = malloc(sizeof(Expr*) * nt);
        for (size_t i = 0; i < nt; i++) t[i] = gr_c(cr[3 * i + j]);
        N[j] = gr_ev(gr_fn(SYM_Plus, t, nt));
        free(t);
    }
    Expr** w = malloc(sizeof(Expr*) * nt);
    Expr** cen = malloc(sizeof(Expr*) * nt);
    for (size_t i = 0; i < nt; i++) {
        w[i] = gr_ev(gr_add(gr_add(gr_mul(gr_c(cr[3 * i]), gr_c(N[0])),
                                   gr_mul(gr_c(cr[3 * i + 1]), gr_c(N[1]))),
                            gr_mul(gr_c(cr[3 * i + 2]), gr_c(N[2]))));
        const Expr* a = gr_at(p, i + 1);
        const Expr* b = gr_at(p, i + 2);
        Expr* c3[3];
        for (size_t j = 0; j < 3; j++)
            c3[j] = gr_mul(gr_fn2(SYM_Rational, expr_new_integer(1), expr_new_integer(3)),
                           gr_add(gr_add(gr_c(gr_at(p0, j)), gr_c(gr_at(a, j))), gr_c(gr_at(b, j))));
        cen[i] = gr_ev(gr_list(c3, 3));
    }
    Expr* out = gr_weighted_mean((const Expr* const*)cen, w, nt, 3);
    for (size_t i = 0; i < nt; i++) { expr_free(w[i]); expr_free(cen[i]); }
    for (size_t i = 0; i < 3 * nt; i++) expr_free(cr[i]);
    for (size_t j = 0; j < 3; j++) expr_free(N[j]);
    free(w); free(cen); free(cr);
    return out;
}

static Expr* gr_centroid_more(const Expr* g) {
    if (g->type == EXPR_FUNCTION && g->data.function.head->type == EXPR_SYMBOL
        && gr_len(g) == 1 && is_ndarray(gr_at(g, 0))) {
        /* packed coordinates: read them as the nested List */
        Expr* u = gr_fn(g->data.function.head->data.symbol.name,
                        (Expr*[]){ gr_unpack(gr_at(g, 0)) }, 1);
        Expr* r = head_is(gr_at(u, 0), SYM_List) ? gr_centroid_more(u) : NULL;
        expr_free(u);
        return r;
    }
    if (head_is(g, SYM_Point)) return gr_centroid_point(g);
    if (head_is(g, SYM_Line)) return gr_centroid_line(g);
    if (head_is(g, SYM_Circle)) return gr_centroid_circle(g);
    if (head_is(g, SYM_Ball) || head_is(g, SYM_Sphere) || head_is(g, SYM_Cuboid))
        return gr_centroid_solid(g);
    return gr_centroid_poly3(g);
}

/* ---- the dispatch ------------------------------------------------------- */

static bool gr_is_ruled_polygon(const Expr* g) {
    return head_is(g, SYM_Polygon) && gr_len(g) == 1 && head_is(gr_at(g, 0), SYM_Rule)
           && gr_len(gr_at(g, 0)) == 2;
}

/* Area/Perimeter/RegionCentroid of a multi-polygon or a polygon with holes. */
static Expr* gr_compound(const char* h, const Expr* g) {
    const Expr* outer = NULL;
    const Expr* parts;
    bool holes = gr_is_ruled_polygon(g);
    if (holes) {
        outer = gr_at(gr_at(g, 0), 0);
        parts = gr_at(gr_at(g, 0), 1);
        size_t d;
        if (!gr_is_point_list(outer, &d)) return NULL;
        if (gr_is_point_list(parts, &d)) parts = NULL;       /* one hole, unwrapped */
        else if (!gr_is_multi(parts)) return NULL;
    } else {
        if (!(head_is(g, SYM_Polygon) && gr_len(g) == 1 && gr_is_multi(gr_at(g, 0)))) return NULL;
        parts = gr_at(g, 0);
    }
    const Expr* hole1 = holes && !parts ? gr_at(gr_at(g, 0), 1) : NULL;
    size_t np = parts ? gr_len(parts) : 1;
    size_t total = np + (holes ? 1 : 0);
    Expr** meas = calloc(total, sizeof(Expr*));
    Expr** cent = calloc(total, sizeof(Expr*));
    bool ok = true;
    bool want_c = strcmp(h, SYM_RegionCentroid) == 0;
    const char* mh = want_c ? SYM_Area : h;
    for (size_t i = 0; i < total && ok; i++) {
        const Expr* pts = holes ? (i == 0 ? outer : (parts ? gr_at(parts, i - 1) : hole1))
                                : gr_at(parts, i);
        meas[i] = gr_call(mh, gr_fn1(SYM_Polygon, gr_c(pts)), NULL);
        if (!meas[i] || !expr_is_numeric_quantity(meas[i])) ok = false;
        if (ok && want_c) {
            cent[i] = gr_call(SYM_RegionCentroid, gr_fn1(SYM_Polygon, gr_c(pts)), NULL);
            if (!cent[i]) ok = false;
        }
    }
    Expr* out = NULL;
    if (ok) {
        /* holes subtract from the area (and from the area-weighted centroid
         * sum), and add to the perimeter */
        bool sub = holes && !strcmp(mh, SYM_Area);
        Expr** w = malloc(sizeof(Expr*) * total);
        Expr** wc = malloc(sizeof(Expr*) * total);
        for (size_t i = 0; i < total; i++) {
            Expr* m = gr_c(meas[i]);
            if (sub && i > 0) m = gr_mul(expr_new_integer(-1), m);
            if (want_c) wc[i] = gr_mul(gr_c(m), gr_c(cent[i]));
            w[i] = m;
        }
        Expr* sum = gr_ev(gr_fn(SYM_Plus, w, total));
        if (want_c) {
            Expr* num = gr_ev(gr_fn(SYM_Plus, wc, total));
            out = gr_ev(gr_fn1(SYM_Simplify, gr_mul(num, gr_inv(sum))));
        } else {
            out = sum;
        }
        free(w); free(wc);
    }
    for (size_t i = 0; i < total; i++) { if (meas[i]) expr_free(meas[i]); if (cent[i]) expr_free(cent[i]); }
    free(meas); free(cent);
    return out;
}

static Expr* gr_member_compound(const Expr* g, const Expr* q) {
    if (gr_is_ruled_polygon(g)) {
        const Expr* outer = gr_at(gr_at(g, 0), 0);
        const Expr* hs = gr_at(gr_at(g, 0), 1);
        size_t d;
        if (!gr_is_point_list(outer, &d)) return NULL;
        Expr* r = gr_call(SYM_RegionMember, gr_fn1(SYM_Polygon, gr_c(outer)), q);
        int in = r ? gr_truth(r) : -1;
        if (in != 1) return in == 0 ? expr_new_symbol(SYM_False) : NULL;
        bool one = gr_is_point_list(hs, &d);
        if (!one && !gr_is_multi(hs)) return NULL;
        size_t nh = one ? 1 : gr_len(hs);
        for (size_t i = 0; i < nh; i++) {
            const Expr* h = one ? hs : gr_at(hs, i);
            Expr* m = gr_call(SYM_RegionMember, gr_fn1(SYM_Polygon, gr_c(h)), q);
            int t = m ? gr_truth(m) : -1;
            if (t < 0) return NULL;
            if (t == 1) {
                /* inside the hole: still a member when on its edge */
                size_t n = gr_len(h);
                Expr** cl = malloc(sizeof(Expr*) * (n + 1));
                for (size_t j = 0; j < n; j++) cl[j] = gr_c(gr_at(h, j));
                cl[n] = gr_c(gr_at(h, 0));
                Expr* ring = gr_list(cl, n + 1);
                free(cl);
                int e = gr_member_line(ring, q);
                expr_free(ring);
                if (e < 0) return NULL;
                if (e == 0) return expr_new_symbol(SYM_False);
            }
        }
        return expr_new_symbol(SYM_True);
    }
    if (head_is(g, SYM_Polygon) && gr_len(g) == 1 && gr_is_multi(gr_at(g, 0))) {
        const Expr* parts = gr_at(g, 0);
        int undecided = 0;
        for (size_t i = 0; i < gr_len(parts); i++) {
            Expr* m = gr_call(SYM_RegionMember, gr_fn1(SYM_Polygon, gr_c(gr_at(parts, i))), q);
            int t = m ? gr_truth(m) : -1;
            if (t == 1) return expr_new_symbol(SYM_True);
            if (t < 0) undecided = 1;
        }
        return undecided ? NULL : expr_new_symbol(SYM_False);
    }
    return NULL;
}

/* The Polygon-with-irrational-coordinates path, Area/Perimeter/Centroid. */
static Expr* gr_symbolic(const char* h, const Expr* g) {
    if (!head_is(g, SYM_Polygon) || gr_len(g) != 1) return NULL;
    const Expr* pts = gr_at(g, 0);
    size_t d;
    if (!gr_is_point_list(pts, &d) || d != 2) return NULL;
    if (gr_all_rational(pts)) return NULL;              /* the GMP path */
    if (gr_has_inexact(pts)) {
        /* Real mixed with Sqrt[2]: numericize, the machine path answers */
        Expr* np = gr_ev(gr_fn1(SYM_N, gr_c(pts)));
        if (expr_eq(np, pts)) { expr_free(np); return NULL; }   /* nothing to numericize */
        return gr_call(h, gr_fn1(SYM_Polygon, np), NULL);
    }
    if (!strcmp(h, SYM_Area)) return gr_sym_area(pts);
    if (!strcmp(h, SYM_Perimeter)) return gr_sym_perimeter(pts);
    if (!strcmp(h, SYM_RegionCentroid)) return gr_sym_centroid(pts);
    return NULL;
}

Expr* geom_region_extra(const char* h, const Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    bool member = !strcmp(h, SYM_RegionMember);
    if (gr_len(res) != (member ? 2u : 1u)) return NULL;
    const Expr* g = gr_at(res, 0);

    Expr* poly = gr_as_polygon(g);
    if (poly) return gr_call(h, poly, member ? gr_at(res, 1) : NULL);

    if (!member) {
        if (!strcmp(h, SYM_RegionCentroid)) {
            Expr* c = gr_centroid_more(g);
            if (c) return c;
        }
        Expr* r = gr_compound(h, g);
        if (r) return r;
        return gr_symbolic(h, g);
    }

    Expr* q = gr_unpack(gr_at(res, 1));
    Expr* out = NULL;
    size_t d;
    if (gr_is_point_list(q, &d) && !head_is(g, SYM_Polygon)) {
        /* a list of points: test each (Polygon already maps its own) */
        size_t k = gr_len(q);
        Expr** items = malloc(sizeof(Expr*) * k);
        for (size_t i = 0; i < k; i++) items[i] = gr_fn2(SYM_RegionMember, gr_c(g), gr_c(gr_at(q, i)));
        out = gr_ev(gr_list(items, k));
        free(items);
        expr_free(q);
        return out;
    }
    if (!gr_is_point(q, 0)) { expr_free(q); return NULL; }
    out = gr_member_compound(g, q);
    if (!out) {
        int t = gr_member_point(g, q);
        if (t >= 0) out = expr_new_symbol(t ? SYM_True : SYM_False);
    }
    if (!out && head_is(g, SYM_Polygon) && gr_len(g) == 1 && gr_exact_points(gr_at(g, 0))
        && gr_len(q) == 2 && !gr_has_inexact(q)
        && !(gr_all_rational(gr_at(g, 0)) && gr_all_rational(q))) {
        int t = gr_sym_member(gr_at(g, 0), q);
        if (t >= 0) out = expr_new_symbol(t ? SYM_True : SYM_False);
    }
    expr_free(q);
    return out;
}
