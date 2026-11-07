/*
 * geometry_curves.c -- ArcLength, RegionMeasure, and the curved regions
 * Annulus / Disk / Circle for the geometry heads (src/geometry.c).
 *
 *   ArcLength[Line[{p1, ..., pn}]]          Sum |p_{i+1} - p_i| (any dimension;
 *                                           a list of point lists sums them)
 *   ArcLength[Circle[c, r]]                 2 Pi r; an ellipse Circle[c, {a, b}]
 *                                           4 max EllipticE[1 - (min/max)^2];
 *                                           an arc Circle[c, r, {t1, t2}] r (t2 - t1)
 *   ArcLength[{x1[t], ..., xn[t]}, {t, a, b}]   Int_a^b |c'(t)| dt
 *   ArcLength[f[t], {t, a, b}]              the graph: Int_a^b Sqrt[1 + f'(t)^2] dt
 *   ArcLength[Point[..]] / ArcLength[Disk[..]]  Undefined (no 1-D measure)
 *
 * The parametric forms integrate the simplified speed |c'(t)|: exactly by
 * Integrate when the range is exact (the closed form Simplify'd), numerically
 * by NIntegrate when a bound is inexact or WorkingPrecision is given -- and
 * also when Integrate finds no closed form over a numeric range.  That last
 * is a deliberate difference: the reference returns the inert Integrate[...]
 * there, but N[] of an inert Integrate is no number in Mathilda, so the inert
 * form would be a dead end.  A symbolic range keeps the inert integral.  A Line whose coordinates are machine
 * numbers sums hypot() over the doubles directly (a contour's thousands of
 * segments never build an expression); exact coordinates give the exact sum of
 * square roots; WorkingPrecision -> p gives N[sum, p].
 *
 *   Annulus[c, {r1, r2}]  (Annulus[] = Annulus[{0, 0}, {1/2, 1}])
 *   Annulus[c, {r1, r2}, {t1, t2}]  the annular sector
 *     Area / RegionMeasure   (t2 - t1)/2 (r2^2 - r1^2), Pi (r2^2 - r1^2) full
 *     Perimeter              2 Pi (r1 + r2) full; (t2 - t1)(r1 + r2) + 2 (r2 - r1)
 *     RegionCentroid         c (full annulus)
 *     RegionMember[A, p]     r1 <= |p - c| <= r2 (closed), and the angle of
 *                            p - c in [t1, t2] for a sector; a symbolic point
 *                            gives the condition (x | y) in Reals && r1^2 <=
 *                            (x - cx)^2 + (y - cy)^2 <= r2^2
 *   Disk[c, r] / Disk[c, {a, b}]: Area / RegionMeasure Pi r^2 / Pi a b.
 *   RegionMeasure of a Polygon is its Area, of a Line or Circle its ArcLength.
 *
 * Numeric radii swapped out of order (Annulus[c, {2, 1}]) measure as the
 * ordered annulus, as in the reference.  None of these is an element-wise
 * numeric map, so no NDArray / Compile surface (documented EXEMPT); the machine
 * Line sum is the one buffer-shaped path and reads packed points directly.
 */

#include "geometry.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "attr.h"
#include "common.h"       /* head_is */
#include "eval.h"
#include "expr.h"
#include "linalg/numarray.h"
#include "ndarray.h"       /* is_ndarray, ndarray_to_nested_list */
#include "numeric.h"
#include "message.h"
#include "print.h"
#include "sym_names.h"
#include "symtab.h"

/* ---- expression helpers (each consumes its arguments) ----------------- */

static Expr* gc_fn(const char* head, Expr** a, size_t n) {
    return expr_new_function(expr_new_symbol(head), a, n);
}
static Expr* gc_fn1(const char* h, Expr* a) { return gc_fn(h, (Expr*[]){ a }, 1); }
static Expr* gc_fn2(const char* h, Expr* a, Expr* b) { return gc_fn(h, (Expr*[]){ a, b }, 2); }
static Expr* gc_add(Expr* a, Expr* b) { return gc_fn2(SYM_Plus, a, b); }
static Expr* gc_mul(Expr* a, Expr* b) { return gc_fn2(SYM_Times, a, b); }
static Expr* gc_neg(Expr* a) { return gc_mul(expr_new_integer(-1), a); }
static Expr* gc_sub(Expr* a, Expr* b) { return gc_add(a, gc_neg(b)); }
static Expr* gc_sq(Expr* a) { return gc_fn2(SYM_Power, a, expr_new_integer(2)); }
static Expr* gc_sqrt(Expr* a) { return gc_fn1("Sqrt", a); }
static Expr* gc_c(const Expr* e) { return expr_copy((Expr*)e); }
static Expr* gc_pi(void) { return expr_new_symbol(SYM_Pi); }
static Expr* gc_ev(Expr* e) { return eval_and_free(e); }

static size_t gc_len(const Expr* e) { return e->data.function.arg_count; }
static const Expr* gc_at(const Expr* e, size_t i) { return e->data.function.args[i]; }

static bool gc_is_inexact_leaf(const Expr* e) {
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* Does `e` contain an inexact number anywhere?  (Called on plain trees only:
 * gc_unpack has already turned any packed / NDArray value into Lists.) */
static bool gc_has_inexact(const Expr* e) {
    if (gc_is_inexact_leaf(e)) return true;
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < gc_len(e); i++) if (gc_has_inexact(gc_at(e, i))) return true;
    return false;
}

/* Does `e` contain an MPFR leaf (an arbitrary-precision computation)? */
static bool gc_has_mpfr(const Expr* e) {
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < gc_len(e); i++) if (gc_has_mpfr(gc_at(e, i))) return true;
    return false;
}

/* ---- options ---------------------------------------------------------- */

/* Strip a trailing WorkingPrecision -> p from the argument list; returns the
 * positional count and sets *wp (borrowed) or leaves it NULL. False on any
 * other option (the call stays unevaluated). */
static bool gc_options(const Expr* res, size_t* npos, const Expr** wp) {
    size_t n = gc_len(res);
    *wp = NULL;
    while (n > 0) {
        const Expr* a = gc_at(res, n - 1);
        if (!(head_is(a, SYM_Rule) || head_is(a, SYM_RuleDelayed)) || gc_len(a) != 2) break;
        const Expr* k = gc_at(a, 0);
        if (k->type == EXPR_SYMBOL && k->data.symbol.name == SYM_WorkingPrecision) {
            *wp = gc_at(a, 1);
            n--;
            continue;
        }
        return false;
    }
    *npos = n;
    return true;
}

static bool gc_wp_is_machine(const Expr* wp) {
    return wp && wp->type == EXPR_SYMBOL && strcmp(wp->data.symbol.name, "MachinePrecision") == 0;
}

/* N[e, wp] for an explicit precision, N[e] for MachinePrecision; consumes e. */
static Expr* gc_apply_wp(Expr* e, const Expr* wp) {
    if (!wp) return e;
    if (gc_wp_is_machine(wp)) return gc_ev(gc_fn1(SYM_N, e));
    return gc_ev(gc_fn2(SYM_N, e, gc_c(wp)));
}

/* ---- Line ------------------------------------------------------------- */

/* Machine sum of segment lengths of one point list; false if the list is not
 * a numeric rank-2 array. */
static bool gc_line_machine(const Expr* pts, double* out) {
    int rows = 0, cols = 0;
    double* buf = NULL;
    if (!na_load_matrix(pts, false, false, &rows, &cols, &buf)) return false;
    double s = 0.0;
    for (int i = 1; i < rows; i++) {
        double ss = 0.0, scale = 0.0;
        /* scaled sum of squares: no overflow for 1e200-sized coordinates */
        for (int j = 0; j < cols; j++) {
            double d = fabs(buf[i * cols + j] - buf[(i - 1) * cols + j]);
            if (d > scale) { ss = ss * (scale / d) * (scale / d) + 1.0; scale = d; }
            else if (d > 0) ss += (d / scale) * (d / scale);
        }
        s += scale * sqrt(ss);
    }
    free(buf);
    *out = s;
    return true;
}

/* Exact (or symbolic, or MPFR) sum of Sqrt[|p_{i+1} - p_i|^2]; NULL when the
 * points are not a list of equal-length coordinate lists. */
static Expr* gc_line_exact(const Expr* pts) {
    size_t n = gc_len(pts);
    if (n == 0) return NULL;
    size_t d = 0;
    for (size_t i = 0; i < n; i++) {
        const Expr* p = gc_at(pts, i);
        if (!head_is(p, SYM_List)) return NULL;
        if (i == 0) d = gc_len(p);
        else if (gc_len(p) != d) return NULL;
        for (size_t j = 0; j < d; j++) if (head_is(gc_at(p, j), SYM_List)) return NULL;
    }
    if (d == 0) return NULL;
    Expr** terms = malloc(sizeof(Expr*) * (n > 1 ? n - 1 : 1));
    for (size_t i = 1; i < n; i++) {
        Expr** sq = malloc(sizeof(Expr*) * d);
        for (size_t j = 0; j < d; j++)
            sq[j] = gc_sq(gc_sub(gc_c(gc_at(gc_at(pts, i), j)), gc_c(gc_at(gc_at(pts, i - 1), j))));
        terms[i - 1] = gc_sqrt(gc_fn(SYM_Plus, sq, d));
        free(sq);
    }
    Expr* sum = (n > 1) ? gc_fn(SYM_Plus, terms, n - 1) : expr_new_integer(0);
    free(terms);
    return gc_ev(sum);
}

/* Is `pts` a list of point LISTS (the multi-line form)? */
static bool gc_is_multiline(const Expr* pts) {
    if (!head_is(pts, SYM_List) || gc_len(pts) == 0) return false;
    const Expr* p0 = gc_at(pts, 0);
    return head_is(p0, SYM_List) && gc_len(p0) > 0 && head_is(gc_at(p0, 0), SYM_List);
}

static Expr* gc_line_length_one(const Expr* pts) {
    if (gc_has_inexact(pts) && !gc_has_mpfr(pts)) {
        double v;
        if (gc_line_machine(pts, &v)) return expr_new_real(v);
    }
    return gc_line_exact(pts);
}

/* A packed / NDArray point list (an EXPR_NDARRAY node, not a List function)
 * as nested Lists, so the exact-vs-machine decision reads the real element
 * heads (an int64 buffer is exact).  Returns a fresh copy either way. */
static Expr* gc_unpack(const Expr* e) {
    if (is_ndarray(e)) {
        Expr* m = ndarray_to_nested_list(e);
        if (m) return m;
    }
    return expr_copy((Expr*)e);
}

static Expr* gc_line_length(const Expr* line) {
    if (gc_len(line) != 1) return NULL;
    Expr* pts = gc_unpack(gc_at(line, 0));
    Expr* out = NULL;
    if (!head_is(pts, SYM_List)) { expr_free(pts); return NULL; }
    if (gc_is_multiline(pts)) {
        /* a list of point lists: the lengths add */
        size_t m = gc_len(pts);
        Expr** parts = malloc(sizeof(Expr*) * m);
        bool ok = true;
        for (size_t i = 0; i < m; i++) {
            parts[i] = gc_line_length_one(gc_at(pts, i));
            if (!parts[i]) {
                for (size_t j = 0; j < i; j++) expr_free(parts[j]);
                ok = false;
                break;
            }
        }
        if (ok) out = gc_ev(gc_fn(SYM_Plus, parts, m));
        free(parts);
    } else {
        out = gc_line_length_one(pts);
    }
    expr_free(pts);
    return out;
}

/* ---- Circle / Disk / Annulus ------------------------------------------ */

/* Radius spec of Circle/Disk: arg 1, default 1. */
static const Expr* gc_radius(const Expr* r) {
    return gc_len(r) >= 2 ? gc_at(r, 1) : NULL;
}

/* 4 b EllipticE[1 - a^2/b^2] with b the larger semi-axis when they are
 * comparable numbers, b = the second otherwise. */
static Expr* gc_ellipse_perimeter(const Expr* a, const Expr* b) {
    const Expr* lo = a; const Expr* hi = b;
    Expr* t = gc_ev(gc_fn2(SYM_Greater, gc_c(a), gc_c(b)));
    if (t->type == EXPR_SYMBOL && t->data.symbol.name == SYM_True) { lo = b; hi = a; }
    expr_free(t);
    Expr* m = gc_sub(expr_new_integer(1), gc_mul(gc_sq(gc_c(lo)), gc_fn2(SYM_Power, gc_sq(gc_c(hi)), expr_new_integer(-1))));
    return gc_ev(gc_mul(gc_mul(expr_new_integer(4), gc_c(hi)), gc_fn1(SYM_EllipticE, m)));
}

static Expr* gc_circle_length(const Expr* c) {
    size_t n = gc_len(c);
    if (n > 3) return NULL;
    const Expr* r = gc_radius(c);
    if (n == 3) {                                       /* arc */
        const Expr* th = gc_at(c, 2);
        if (!head_is(th, SYM_List) || gc_len(th) != 2 || head_is(r, SYM_List)) return NULL;
        return gc_ev(gc_mul(gc_c(r), gc_sub(gc_c(gc_at(th, 1)), gc_c(gc_at(th, 0)))));
    }
    if (!r) return gc_ev(gc_mul(expr_new_integer(2), gc_pi()));
    if (head_is(r, SYM_List)) {
        if (gc_len(r) != 2) return NULL;
        return gc_ellipse_perimeter(gc_at(r, 0), gc_at(r, 1));
    }
    return gc_ev(gc_mul(gc_mul(expr_new_integer(2), gc_pi()), gc_c(r)));
}

static Expr* gc_disk_area(const Expr* d) {
    size_t n = gc_len(d);
    if (n > 2) return NULL;                             /* sectors: not handled */
    const Expr* r = gc_radius(d);
    if (!r) return gc_pi();
    if (head_is(r, SYM_List)) {
        if (gc_len(r) != 2) return NULL;
        return gc_ev(gc_mul(gc_pi(), gc_mul(gc_c(gc_at(r, 0)), gc_c(gc_at(r, 1)))));
    }
    return gc_ev(gc_mul(gc_pi(), gc_sq(gc_c(r))));
}

/* Unpack Annulus[c, {r1, r2}(, {t1, t2})]: borrowed parts; r1 <= r2 when both
 * are comparable numbers. */
typedef struct { const Expr *c, *r1, *r2, *t1, *t2; } GcAnnulus;

static bool gc_annulus(const Expr* a, GcAnnulus* out) {
    if (!head_is(a, SYM_Annulus)) return false;
    size_t n = gc_len(a);
    if (n < 2 || n > 3) return false;
    const Expr* c = gc_at(a, 0);
    const Expr* rr = gc_at(a, 1);
    if (!head_is(c, SYM_List) || gc_len(c) != 2) return false;
    if (!head_is(rr, SYM_List) || gc_len(rr) != 2) return false;
    out->c = c;
    out->r1 = gc_at(rr, 0);
    out->r2 = gc_at(rr, 1);
    out->t1 = out->t2 = NULL;
    if (n == 3) {
        const Expr* th = gc_at(a, 2);
        if (!head_is(th, SYM_List) || gc_len(th) != 2) return false;
        out->t1 = gc_at(th, 0);
        out->t2 = gc_at(th, 1);
    }
    Expr* t = gc_ev(gc_fn2(SYM_Greater, gc_c(out->r1), gc_c(out->r2)));
    if (t->type == EXPR_SYMBOL && t->data.symbol.name == SYM_True) {
        const Expr* s = out->r1; out->r1 = out->r2; out->r2 = s;
    }
    expr_free(t);
    return true;
}

static Expr* gc_annulus_area(const GcAnnulus* A) {
    Expr* d = gc_sub(gc_sq(gc_c(A->r2)), gc_sq(gc_c(A->r1)));
    if (!A->t1) return gc_ev(gc_mul(gc_pi(), d));
    Expr* half = gc_fn2(SYM_Rational, expr_new_integer(1), expr_new_integer(2));
    return gc_ev(gc_mul(gc_mul(half, gc_sub(gc_c(A->t2), gc_c(A->t1))), d));
}

static Expr* gc_annulus_perimeter(const GcAnnulus* A) {
    Expr* s = gc_add(gc_c(A->r1), gc_c(A->r2));
    if (!A->t1) return gc_ev(gc_mul(gc_mul(expr_new_integer(2), gc_pi()), s));
    Expr* arcs = gc_mul(gc_sub(gc_c(A->t2), gc_c(A->t1)), s);
    Expr* radial = gc_mul(expr_new_integer(2), gc_sub(gc_c(A->r2), gc_c(A->r1)));
    return gc_ev(gc_add(arcs, radial));
}

/* The region's 2-D measure, or NULL. */
static Expr* gc_area_of(const Expr* region) {
    GcAnnulus A;
    if (gc_annulus(region, &A)) return gc_annulus_area(&A);
    if (head_is(region, SYM_Disk)) return gc_disk_area(region);
    return NULL;
}

/* ---- hooks called from src/geometry.c --------------------------------- */

Expr* geom_curved_area(const Expr* region) { return gc_area_of(region); }

Expr* geom_curved_perimeter(const Expr* region) {
    GcAnnulus A;
    if (gc_annulus(region, &A)) return gc_annulus_perimeter(&A);
    if (head_is(region, SYM_Disk) && gc_len(region) <= 2) {
        /* the boundary circle (or ellipse) */
        Expr* circ = expr_copy((Expr*)region);
        expr_free(circ->data.function.head);
        circ->data.function.head = expr_new_symbol(SYM_Circle);
        Expr* r = gc_circle_length(circ);
        expr_free(circ);
        return r;
    }
    return NULL;
}

Expr* geom_curved_centroid(const Expr* region) {
    GcAnnulus A;
    if (gc_annulus(region, &A) && !A.t1) return gc_c(A.c);
    if (head_is(region, SYM_Disk) && gc_len(region) <= 2) {
        if (gc_len(region) == 0) return gc_ev(gc_fn2(SYM_List, expr_new_integer(0), expr_new_integer(0)));
        return gc_c(gc_at(region, 0));
    }
    return NULL;
}

/* RegionMember[Annulus[..], p]: True/False for a numeric point, the condition
 * for a symbolic one, NULL otherwise. */
static Expr* gc_member_point(const GcAnnulus* A, const Expr* p);

Expr* geom_curved_member(const Expr* region, const Expr* p) {
    GcAnnulus A;
    if (!gc_annulus(region, &A)) return NULL;
    Expr* q = gc_unpack(p);                /* a packed point is an NDArray node */
    Expr* r = gc_member_point(&A, q);
    expr_free(q);
    return r;
}

static Expr* gc_member_point(const GcAnnulus* Ap, const Expr* p) {
    GcAnnulus A = *Ap;
    if (!head_is(p, SYM_List) || gc_len(p) != 2) return NULL;
    Expr* dx = gc_sub(gc_c(gc_at(p, 0)), gc_c(gc_at(A.c, 0)));
    Expr* dy = gc_sub(gc_c(gc_at(p, 1)), gc_c(gc_at(A.c, 1)));
    Expr* rho2 = gc_ev(gc_add(gc_sq(dx), gc_sq(dy)));
    Expr* radial = gc_ev(gc_fn(SYM_LessEqual,
        (Expr*[]){ gc_ev(gc_sq(gc_c(A.r1))), gc_c(rho2), gc_ev(gc_sq(gc_c(A.r2))) }, 3));
    bool numeric = expr_is_numeric_quantity(gc_at(p, 0)) && expr_is_numeric_quantity(gc_at(p, 1));
    if (!numeric) {
        expr_free(rho2);
        if (A.t1) { expr_free(radial); return NULL; }
        /* (x | y) in Reals && r1^2 <= rho^2 <= r2^2 */
        Expr* alts = gc_fn2(SYM_Alternatives, gc_c(gc_at(p, 0)), gc_c(gc_at(p, 1)));
        Expr* elem = gc_fn2(SYM_Element, alts, expr_new_symbol(SYM_Reals));
        return gc_ev(gc_fn2(SYM_And, elem, radial));
    }
    bool in = radial->type == EXPR_SYMBOL && radial->data.symbol.name == SYM_True;
    bool decided = radial->type == EXPR_SYMBOL &&
                   (radial->data.symbol.name == SYM_True || radial->data.symbol.name == SYM_False);
    expr_free(radial);
    if (!decided) { expr_free(rho2); return NULL; }
    if (in && A.t1) {
        /* angle of p - c, taken into [t1, t1 + 2 Pi) */
        double t1, t2, px, py, cx, cy, im;
        bool ok = na_read_scalar(A.t1, &t1, &im) && na_read_scalar(A.t2, &t2, &im);
        Expr* t1n = gc_ev(gc_fn1(SYM_N, gc_c(A.t1)));
        Expr* t2n = gc_ev(gc_fn1(SYM_N, gc_c(A.t2)));
        Expr* pxn = gc_ev(gc_fn1(SYM_N, gc_c(gc_at(p, 0))));
        Expr* pyn = gc_ev(gc_fn1(SYM_N, gc_c(gc_at(p, 1))));
        Expr* cxn = gc_ev(gc_fn1(SYM_N, gc_c(gc_at(A.c, 0))));
        Expr* cyn = gc_ev(gc_fn1(SYM_N, gc_c(gc_at(A.c, 1))));
        ok = na_read_scalar(t1n, &t1, &im) && na_read_scalar(t2n, &t2, &im) &&
             na_read_scalar(pxn, &px, &im) && na_read_scalar(pyn, &py, &im) &&
             na_read_scalar(cxn, &cx, &im) && na_read_scalar(cyn, &cy, &im);
        expr_free(t1n); expr_free(t2n); expr_free(pxn); expr_free(pyn);
        expr_free(cxn); expr_free(cyn);
        if (!ok) { expr_free(rho2); return NULL; }
        if (px != cx || py != cy) {
            double ang = atan2(py - cy, px - cx);
            double two_pi = 2.0 * 3.14159265358979323846;
            while (ang < t1) ang += two_pi;
            while (ang >= t1 + two_pi) ang -= two_pi;
            in = ang <= t2 + 1e-15 * fabs(t2);
        }
    }
    expr_free(rho2);
    return expr_new_symbol(in ? SYM_True : SYM_False);
}

/* ---- ArcLength --------------------------------------------------------- */

static bool gc_contains_head(const Expr* e, const char* h) {
    if (e->type != EXPR_FUNCTION) return false;
    if (head_is(e, h)) return true;
    if (gc_contains_head(e->data.function.head, h)) return true;
    for (size_t i = 0; i < gc_len(e); i++) if (gc_contains_head(gc_at(e, i), h)) return true;
    return false;
}

/* The speed |c'(t)| of a curve c (a List of components, or one function whose
 * graph is meant), simplified; NULL on failure. */
static Expr* gc_speed(const Expr* curve, const Expr* t) {
    Expr* sq;
    if (head_is(curve, SYM_List)) {
        size_t n = gc_len(curve);
        if (n == 0) return NULL;
        Expr** parts = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++)
            parts[i] = gc_sq(gc_fn2(SYM_D, gc_c(gc_at(curve, i)), gc_c(t)));
        sq = gc_fn(SYM_Plus, parts, n);
        free(parts);
    } else {
        sq = gc_add(expr_new_integer(1), gc_sq(gc_fn2(SYM_D, gc_c(curve), gc_c(t))));
    }
    /* Simplify the squared speed first (Sin^2 + Cos^2 -> 1), then the root:
     * Sqrt[Simplify[...]] keeps Sqrt[1 + 4 t^2] rather than PowerExpanding. */
    Expr* s = gc_ev(gc_fn1(SYM_Simplify, gc_ev(sq)));
    return gc_ev(gc_fn1(SYM_Simplify, gc_sqrt(s)));
}

static Expr* gc_arclength_param(const Expr* curve, const Expr* range, const Expr* wp) {
    if (!head_is(range, SYM_List) || gc_len(range) != 3) return NULL;
    const Expr* t = gc_at(range, 0);
    if (t->type != EXPR_SYMBOL) return NULL;
    Expr* speed = gc_speed(curve, t);
    if (!speed) return NULL;
    bool numeric = wp != NULL || gc_has_inexact(gc_at(range, 1)) || gc_has_inexact(gc_at(range, 2));
    if (speed->type == EXPR_INTEGER && speed->data.integer == 0) {
        expr_free(speed);
        return expr_new_integer(0);
    }
    if (!numeric) {
        /* Exactly first, quietly: when Integrate finds no closed form its
         * complaint is noise, because the numeric value follows. */
        Expr* ex = gc_ev(gc_fn1("Quiet", gc_fn2(SYM_Integrate, gc_c(speed), gc_c(range))));
        if (!gc_contains_head(ex, SYM_Integrate)) {
            expr_free(speed);
            return gc_ev(gc_fn1(SYM_Simplify, ex));
        }
        expr_free(ex);
        /* No closed form: the numeric value when the range is numeric (the
         * reference keeps the inert Integrate here, but N[] of an inert
         * Integrate is not a number in Mathilda, so that would be a dead end);
         * a symbolic range keeps the inert integral. */
        if (!expr_is_numeric_quantity(gc_at(range, 1)) || !expr_is_numeric_quantity(gc_at(range, 2)))
            return gc_ev(gc_fn2(SYM_Integrate, speed, gc_c(range)));
    }
    Expr* r = gc_fn(SYM_NIntegrate,
        (Expr*[]){ speed, gc_c(range),
                   gc_fn2(SYM_Rule, expr_new_symbol(SYM_WorkingPrecision),
                          wp ? gc_c(wp) : expr_new_symbol("MachinePrecision")) }, 3);
    return gc_ev(r);
}

Expr* builtin_arclength(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t npos;
    const Expr* wp;
    if (!gc_options(res, &npos, &wp)) return NULL;
    if (npos == 1) {
        const Expr* g = gc_at(res, 0);
        Expr* r = NULL;
        if (head_is(g, SYM_Line)) r = gc_line_length(g);
        else if (head_is(g, SYM_Circle)) r = gc_circle_length(g);
        else if (head_is(g, SYM_Point) || head_is(g, SYM_Disk) || head_is(g, SYM_Annulus))
            return expr_new_symbol(SYM_Undefined);
        return r ? gc_apply_wp(r, wp) : NULL;
    }
    if (npos == 2) return gc_arclength_param(gc_at(res, 0), gc_at(res, 1), wp);
    return NULL;
}

/* ---- RegionMeasure ------------------------------------------------------ */

Expr* builtin_regionmeasure(Expr* res) {
    if (res->type != EXPR_FUNCTION || gc_len(res) != 1) return NULL;
    const Expr* g = gc_at(res, 0);
    if (head_is(g, SYM_Polygon)) {
        Expr* call = gc_fn1(SYM_Area, gc_c(g));
        Expr* r = gc_ev(call);
        if (head_is(r, SYM_Area)) { expr_free(r); return NULL; }
        return r;
    }
    if (head_is(g, SYM_Line) || head_is(g, SYM_Circle)) {
        Expr* r = builtin_arclength(res);
        return r;
    }
    return gc_area_of(g);
}

/* ---- Annulus ------------------------------------------------------------ */

/* Annulus[] = Annulus[{0, 0}, {1/2, 1}]; Annulus[c] = Annulus[c, {1/2, 1}]. */
Expr* builtin_annulus(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = gc_len(res);
    if (n > 1) return NULL;
    Expr* c = n == 1 ? gc_c(gc_at(res, 0))
                     : gc_fn2(SYM_List, expr_new_integer(0), expr_new_integer(0));
    Expr* r = gc_fn2(SYM_List, gc_fn2(SYM_Rational, expr_new_integer(1), expr_new_integer(2)),
                     expr_new_integer(1));
    return gc_ev(gc_fn2(SYM_Annulus, c, r));
}

void geometry_curves_init(void) {
    symtab_add_builtin("ArcLength", builtin_arclength);
    symtab_add_builtin("RegionMeasure", builtin_regionmeasure);
    symtab_add_builtin("Annulus", builtin_annulus);
    symtab_get_def("ArcLength")->attributes |= ATTR_PROTECTED;
    symtab_get_def("RegionMeasure")->attributes |= ATTR_PROTECTED;
    symtab_get_def("Annulus")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ArcLength",
        "ArcLength[Line[{p1, ..., pn}]] gives the length of the polyline (any dimension).\n"
        "ArcLength[Circle[c, r]] gives 2 Pi r (an ellipse or an arc likewise).\n"
        "ArcLength[{x1, ..., xn}, {t, a, b}] gives the length of the parametrised curve.\n"
        "ArcLength[f, {t, a, b}] gives the length of the graph of f.\n"
        "Exact ranges integrate exactly (an inert Integrate when no closed form is found);\n"
        "an inexact bound or WorkingPrecision -> p integrates numerically. Point and Disk\n"
        "have Undefined arc length.");
    symtab_set_docstring("RegionMeasure",
        "RegionMeasure[reg] gives the measure of a region: the area of a Polygon,\n"
        "Disk or Annulus, the length of a Line or Circle.");
    symtab_set_docstring("Annulus",
        "Annulus[{x, y}, {rin, rout}] represents the closed annulus rin <= |p - {x, y}| <= rout.\n"
        "Annulus[{x, y}, {rin, rout}, {t1, t2}] represents the annular sector t1 <= angle <= t2.\n"
        "Annulus[] is Annulus[{0, 0}, {1/2, 1}]. Area, RegionMeasure, Perimeter,\n"
        "RegionCentroid and RegionMember accept it.");
}
