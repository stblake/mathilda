/*
 * geometry.c -- core computational-geometry builtins (GEO-1).
 *
 * Area / Perimeter / RegionCentroid / RegionMember over 2D Polygon, plus
 * ConvexHullRegion over a 2D point set. Semantics follow Wolfram Language 12
 * (reference outputs verified against a live kernel; see
 * docs/spec/builtins/geometry.md and thoughts/shared/tickets/GEO-1/research.md).
 *
 * Two computation paths per head, chosen once per call:
 *
 *   EXACT   -- every coordinate is Integer / BigInt / Rational: all arithmetic
 *              in GMP mpq_t, results canonicalized via make_rational_mpz, so
 *              Area[Polygon[{{0,0},{1,0},{1/2,1/2}}]] is exactly 1/4 and
 *              Perimeter builds symbolic Plus[Sqrt[...],...] trees the
 *              evaluator simplifies (2 + Sqrt[2], 3/2 + 1/2 Sqrt[5], ...).
 *   MACHINE -- any Real/MPFR coordinate switches the whole computation to
 *              doubles (WL's contagion semantics). A visible NDArray points
 *              argument is read via na_load_matrix and is always machine.
 *
 * Anything out of scope -- symbolic coordinates, non-2D data, Polygon with
 *  holes, self-intersecting semantics -- declines (returns NULL, never freeing
 * `res`), leaving the expression unevaluated per the builtin contract.
 *
 * Deliberate, documented deviations from WL (spec page carries the list):
 *   - RegionCentroid declines on zero-area (degenerate) polygons instead of
 *     returning the lower-dimensional measure centroid.
 *   - ConvexHullRegion returns bare Polygon[{verts}] without WL's cell-spec
 *     second argument.
 *   - Machine-path boundary membership and the zero-area gate use exact
 *     IEEE comparisons (no epsilon).
 *
 * Simple (non-self-intersecting) polygons only: Area/RegionCentroid on a
 * self-intersecting vertex list follow shoelace/even-odd semantics, which
 * differ from WL's enclosed-region model. Not detected at runtime.
 */

#include "geometry.h"

#include <gmp.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "arithmetic.h"   /* make_rational_mpz */
#include "attr.h"
#include "common.h"       /* head_is */
#include "eval.h"         /* eval_and_free */
#include "expr.h"
#include "ndarray.h"      /* is_ndarray */
#include "numarray.h"     /* na_read_scalar, na_load_vector, na_load_matrix */
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------------- */
/* Point-list loading                                                        */
/* ------------------------------------------------------------------------- */

typedef struct {
    size_t n;
    int    exact;   /* 1: qx/qy valid (and dx/dy mirror them); 0: dx/dy only */
    int    d_finite; /* 1 iff every dx/dy is finite; an exact coordinate too
                      * large for a double (10^400) mirrors to inf, and any
                      * machine-path arithmetic on it yields NaN cross products
                      * that read as "collinear". Machine paths must decline
                      * rather than answer from those. */
    mpq_t* qx;
    mpq_t* qy;
    double* dx;     /* always filled, so a machine consumer never re-reads */
    double* dy;
} GeomPoints;

static void geom_points_clear(GeomPoints* p) {
    if (p->exact && p->qx) {
        for (size_t i = 0; i < p->n; i++) { mpq_clear(p->qx[i]); mpq_clear(p->qy[i]); }
    }
    free(p->qx); free(p->qy); free(p->dx); free(p->dy);
    memset(p, 0, sizeof *p);
}

/* Exact leaf: Integer, BigInt, or Rational[int-like, int-like]. */
static int geom_leaf_is_exact(const Expr* e) {
    if (expr_is_integer_like(e)) return 1;
    if (head_is(e, SYM_Rational) && e->data.function.arg_count == 2 &&
        expr_is_integer_like(e->data.function.args[0]) &&
        expr_is_integer_like(e->data.function.args[1]))
        return 1;
    return 0;
}

/* Precondition: geom_leaf_is_exact(e). `out` must be initialized.
 *
 * expr_to_mpz INITIALIZES its output (mpz_init_set_*, src/expr.c), so it must
 * NEVER be handed mpq_numref/mpq_denref of a live mpq_t: that overwrites the
 * limb pointer mpq_init already allocated and orphans it. Go through a
 * temporary and copy in. (Found by `leaks --atExit`: one 16-byte block per
 * call with a rational coordinate -- the denominator, the component mpq_init
 * really does pre-allocate.) */
static void geom_leaf_to_mpq(const Expr* e, mpq_t out) {
    mpz_t t;
    if (expr_is_integer_like(e)) {
        expr_to_mpz(e, t);              /* initializes t */
        mpz_set(mpq_numref(out), t);
        mpz_clear(t);
        mpz_set_ui(mpq_denref(out), 1);
        return;
    }
    expr_to_mpz(e->data.function.args[0], t);
    mpz_set(mpq_numref(out), t);
    mpz_clear(t);
    expr_to_mpz(e->data.function.args[1], t);
    mpz_set(mpq_denref(out), t);
    mpz_clear(t);
    mpq_canonicalize(out);
}

/* Machine-numeric leaf (Integer/BigInt/Rational/Real/MPFR, imaginary part 0). */
static int geom_leaf_to_double(const Expr* e, double* out) {
    double re = 0.0, im = 0.0;
    if (!na_read_scalar(e, &re, &im)) return 0;
    if (im != 0.0) return 0;
    if (!isfinite(re)) return 0;   /* 10^400 reads as inf; refuse, never compute on it */
    *out = re;
    return 1;
}

static int geom_alloc(GeomPoints* p, size_t n, int exact) {
    memset(p, 0, sizeof *p);
    p->n = n;
    p->exact = exact;
    p->d_finite = 1;
    p->dx = (double*)malloc(n * sizeof(double));
    p->dy = (double*)malloc(n * sizeof(double));
    if (!p->dx || !p->dy) { geom_points_clear(p); return 0; }
    if (exact) {
        p->qx = (mpq_t*)malloc(n * sizeof(mpq_t));
        p->qy = (mpq_t*)malloc(n * sizeof(mpq_t));
        if (!p->qx || !p->qy) { p->exact = 0; geom_points_clear(p); return 0; }
        for (size_t i = 0; i < n; i++) { mpq_init(p->qx[i]); mpq_init(p->qy[i]); }
    }
    return 1;
}

/* Every mirrored double finite? Sets and returns p->d_finite. */
static int geom_check_finite(GeomPoints* p) {
    for (size_t i = 0; i < p->n; i++) {
        if (!isfinite(p->dx[i]) || !isfinite(p->dy[i])) { p->d_finite = 0; return 0; }
    }
    p->d_finite = 1;
    return 1;
}

/* Load `pts` -- a List of 2D points, or a visible NDArray (rank 2, 2 columns)
 * -- into `out`. Returns 1 on success (caller must geom_points_clear), 0 to
 * decline. A packed List reaching us presents as EXPR_NDARRAY too and takes
 * the same machine path. */
static int geom_read_points(const Expr* pts, GeomPoints* out) {
    memset(out, 0, sizeof *out);

    if (is_ndarray(pts)) {
        /* An INT64-typed ndarray holds exact machine integers, so it must give
         * the same exact answer as the same points written as a nested List --
         * SPEC.md §9 requires the representations to agree. Reading it through
         * the double path would silently downgrade Area[...] from 4 to 4.0. */
        const NDArrayData* nd = &pts->data.ndarray;
        if (nd->rank == 2 && nd->dims[1] == 2 && nd->dims[0] >= 1 &&
            nd->dtype == NDT_INT64) {
            size_t rows = (size_t)nd->dims[0];
            if (!geom_alloc(out, rows, 1)) return 0;
            for (size_t i = 0; i < rows; i++) {
                mpz_set_si(mpq_numref(out->qx[i]), (long)ndt_get_i(nd->data, 2 * i, nd->dtype));
                mpz_set_ui(mpq_denref(out->qx[i]), 1);
                mpz_set_si(mpq_numref(out->qy[i]), (long)ndt_get_i(nd->data, 2 * i + 1, nd->dtype));
                mpz_set_ui(mpq_denref(out->qy[i]), 1);
                out->dx[i] = mpq_get_d(out->qx[i]);
                out->dy[i] = mpq_get_d(out->qy[i]);
            }
            return 1;
        }
        int rows = 0, cols = 0;
        double* buf = NULL;
        if (!na_load_matrix(pts, false, false, &rows, &cols, &buf)) return 0;
        if (cols != 2 || rows < 1) { free(buf); return 0; }
        if (!geom_alloc(out, (size_t)rows, 0)) { free(buf); return 0; }
        for (int i = 0; i < rows; i++) {
            out->dx[i] = buf[2 * i];
            out->dy[i] = buf[2 * i + 1];
        }
        free(buf);
        if (!geom_check_finite(out)) { geom_points_clear(out); return 0; }
        return 1;
    }

    if (!head_is(pts, SYM_List)) return 0;
    size_t n = pts->data.function.arg_count;
    if (n < 1) return 0;

    /* Classification pass: exact until proven machine; decline on symbolic. */
    int exact = 1;
    for (size_t i = 0; i < n; i++) {
        const Expr* pt = pts->data.function.args[i];
        if (!head_is(pt, SYM_List) || pt->data.function.arg_count != 2) return 0;
        for (size_t c = 0; c < 2; c++) {
            const Expr* leaf = pt->data.function.args[c];
            double d;
            if (geom_leaf_is_exact(leaf)) continue;
            if (geom_leaf_to_double(leaf, &d)) { exact = 0; continue; }
            return 0;  /* symbolic / complex / non-numeric: decline */
        }
    }

    if (!geom_alloc(out, n, exact)) return 0;
    for (size_t i = 0; i < n; i++) {
        const Expr* pt = pts->data.function.args[i];
        const Expr* ex = pt->data.function.args[0];
        const Expr* ey = pt->data.function.args[1];
        if (exact) {
            geom_leaf_to_mpq(ex, out->qx[i]);
            geom_leaf_to_mpq(ey, out->qy[i]);
            out->dx[i] = mpq_get_d(out->qx[i]);
            out->dy[i] = mpq_get_d(out->qy[i]);
        } else {
            if (!geom_leaf_to_double(ex, &out->dx[i]) ||
                !geom_leaf_to_double(ey, &out->dy[i])) {
                geom_points_clear(out);
                return 0;
            }
        }
    }
    /* A machine-path point set with a non-finite coordinate cannot be computed
     * on: decline. For an EXACT set the mirror may overflow to inf while the
     * exact arithmetic is still perfectly fine, so only record the fact. */
    geom_check_finite(out);
    if (!out->exact && !out->d_finite) { geom_points_clear(out); return 0; }
    return 1;
}

/* Collapse consecutive duplicate vertices, including the wrap-around pair, so
 * `n` is the count of DISTINCT vertices the docstring and spec promise. This
 * covers the explicitly-closed form ({p1,...,pk,p1}, which WL accepts) and also
 * repeats in the middle: Area[Polygon[{{0,0},{0,0},{1,1}}]] has two distinct
 * vertices and must be Undefined, not a confident 0. */
static int geom_same_point(const GeomPoints* p, size_t i, size_t j) {
    if (p->exact)
        return mpq_equal(p->qx[i], p->qx[j]) && mpq_equal(p->qy[i], p->qy[j]);
    return (p->dx[i] == p->dx[j]) && (p->dy[i] == p->dy[j]);
}

static void geom_dedup_consecutive(GeomPoints* p) {
    if (p->n < 2) return;
    size_t w = 1;
    for (size_t i = 1; i < p->n; i++) {
        if (geom_same_point(p, w - 1, i)) continue;   /* drop the repeat */
        if (w != i) {
            if (p->exact) {
                mpq_set(p->qx[w], p->qx[i]);
                mpq_set(p->qy[w], p->qy[i]);
            }
            p->dx[w] = p->dx[i];
            p->dy[w] = p->dy[i];
        }
        w++;
    }
    /* wrap-around: last equal to first is a repeat too */
    if (w > 1 && geom_same_point(p, w - 1, 0)) w--;
    if (p->exact)
        for (size_t k = w; k < p->n; k++) { mpq_clear(p->qx[k]); mpq_clear(p->qy[k]); }
    p->n = w;
}

/* ------------------------------------------------------------------------- */
/* Small shared helpers                                                      */
/* ------------------------------------------------------------------------- */

/* Canonical exact scalar from an mpq (Integer when the denominator is 1). */
static Expr* geom_mpq_to_expr(const mpq_t q) {
    return make_rational_mpz(mpq_numref(q), mpq_denref(q));
}

static Expr* geom_list2(Expr* a, Expr* b) {
    Expr** args = (Expr**)malloc(2 * sizeof(Expr*));
    if (!args) { expr_free(a); expr_free(b); return NULL; }
    args[0] = a; args[1] = b;
    /* expr_new_function COPIES the args array (src/expr.c:257) rather than
     * taking ownership, so the caller's buffer must be freed here -- the
     * repo-wide convention (see src/list/join.c:50-51). */
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), args, 2);
    free(args);
    return out;
}

/* sign of (b - a) x (c - a), exact. */
static int geom_cross_sign_q(const mpq_t ax, const mpq_t ay,
                             const mpq_t bx, const mpq_t by,
                             const mpq_t cx, const mpq_t cy) {
    mpq_t ux, uy, vx, vy, t1, t2;
    mpq_inits(ux, uy, vx, vy, t1, t2, (mpq_ptr)NULL);
    mpq_sub(ux, bx, ax); mpq_sub(uy, by, ay);
    mpq_sub(vx, cx, ax); mpq_sub(vy, cy, ay);
    mpq_mul(t1, ux, vy); mpq_mul(t2, uy, vx);
    int s = mpq_cmp(t1, t2);
    mpq_clears(ux, uy, vx, vy, t1, t2, (mpq_ptr)NULL);
    return (s > 0) - (s < 0);
}

static int geom_cross_sign_d(double ax, double ay, double bx, double by,
                             double cx, double cy) {
    double cr = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
    return (cr > 0.0) - (cr < 0.0);
}

/* Signed area * 2 (shoelace sum), exact. `out` must be initialized. */
static void geom_signed_area2_q(const GeomPoints* p, mpq_t out) {
    mpq_t t1, t2;
    mpq_inits(t1, t2, (mpq_ptr)NULL);
    mpq_set_ui(out, 0, 1);
    for (size_t i = 0; i < p->n; i++) {
        size_t j = (i + 1) % p->n;
        mpq_mul(t1, p->qx[i], p->qy[j]);
        mpq_mul(t2, p->qx[j], p->qy[i]);
        mpq_sub(t1, t1, t2);
        mpq_add(out, out, t1);
    }
    mpq_clears(t1, t2, (mpq_ptr)NULL);
}

static double geom_signed_area2_d(const GeomPoints* p) {
    double s = 0.0;
    for (size_t i = 0; i < p->n; i++) {
        size_t j = (i + 1) % p->n;
        s += p->dx[i] * p->dy[j] - p->dx[j] * p->dy[i];
    }
    return s;
}

/* Unwrap Area/Perimeter/RegionCentroid's single argument: Polygon[pts] with
 * exactly one part. Returns the points expression or NULL. */
static const Expr* geom_polygon_points(const Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    const Expr* poly = res->data.function.args[0];
    if (!head_is(poly, SYM_Polygon) || poly->data.function.arg_count != 1) return NULL;
    return poly->data.function.args[0];
}

/* ------------------------------------------------------------------------- */
/* Area                                                                      */
/* ------------------------------------------------------------------------- */

Expr* builtin_area(Expr* res) {
    const Expr* pts = geom_polygon_points(res);
    if (!pts) return NULL;

    GeomPoints p;
    if (!geom_read_points(pts, &p)) return NULL;
    geom_dedup_consecutive(&p);

    if (p.n < 3) {  /* degenerate polygon: WL gives Undefined */
        geom_points_clear(&p);
        return expr_new_symbol(SYM_Undefined);
    }

    Expr* out;
    if (p.exact) {
        mpq_t a2;
        mpq_init(a2);
        geom_signed_area2_q(&p, a2);
        mpq_abs(a2, a2);
        mpq_div_2exp(a2, a2, 1);
        out = geom_mpq_to_expr(a2);
        mpq_clear(a2);
    } else {
        out = expr_new_real(fabs(geom_signed_area2_d(&p)) / 2.0);
    }
    geom_points_clear(&p);
    return out;
}

/* ------------------------------------------------------------------------- */
/* Perimeter                                                                 */
/* ------------------------------------------------------------------------- */

Expr* builtin_perimeter(Expr* res) {
    const Expr* pts = geom_polygon_points(res);
    if (!pts) return NULL;

    GeomPoints p;
    if (!geom_read_points(pts, &p)) return NULL;
    geom_dedup_consecutive(&p);

    if (p.n < 3) {
        geom_points_clear(&p);
        return expr_new_symbol(SYM_Undefined);
    }

    Expr* out;
    if (p.exact) {
        /* Sum of Sqrt[dx^2 + dy^2] terms, handed to the evaluator so radical
         * canonicalization (Sqrt[8] -> 2 Sqrt[2], Sqrt[5/4] -> 1/2 Sqrt[5])
         * produces the WL-shaped exact result. */
        Expr** terms = (Expr**)malloc(p.n * sizeof(Expr*));
        if (!terms) { geom_points_clear(&p); return NULL; }
        mpq_t dx, dy, r;
        mpq_inits(dx, dy, r, (mpq_ptr)NULL);
        for (size_t i = 0; i < p.n; i++) {
            size_t j = (i + 1) % p.n;
            mpq_sub(dx, p.qx[j], p.qx[i]);
            mpq_sub(dy, p.qy[j], p.qy[i]);
            mpq_mul(dx, dx, dx);
            mpq_mul(dy, dy, dy);
            mpq_add(r, dx, dy);
            Expr** sq = (Expr**)malloc(sizeof(Expr*));
            if (!sq) {  /* OOM: free the terms built so far and decline */
                for (size_t k = 0; k < i; k++) expr_free(terms[k]);
                free(terms);
                mpq_clears(dx, dy, r, (mpq_ptr)NULL);
                geom_points_clear(&p);
                return NULL;
            }
            sq[0] = geom_mpq_to_expr(r);
            terms[i] = expr_new_function(expr_new_symbol(SYM_Sqrt), sq, 1);
            free(sq);
        }
        mpq_clears(dx, dy, r, (mpq_ptr)NULL);
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, p.n);
        free(terms);
        out = eval_and_free(sum);
    } else {
        double s = 0.0;
        for (size_t i = 0; i < p.n; i++) {
            size_t j = (i + 1) % p.n;
            s += hypot(p.dx[j] - p.dx[i], p.dy[j] - p.dy[i]);
        }
        out = expr_new_real(s);
    }
    geom_points_clear(&p);
    return out;
}

/* ------------------------------------------------------------------------- */
/* RegionCentroid                                                            */
/* ------------------------------------------------------------------------- */

Expr* builtin_region_centroid(Expr* res) {
    const Expr* pts = geom_polygon_points(res);
    if (!pts) return NULL;

    GeomPoints p;
    if (!geom_read_points(pts, &p)) return NULL;
    geom_dedup_consecutive(&p);

    if (p.n < 3) { geom_points_clear(&p); return NULL; }

    Expr* out = NULL;
    if (p.exact) {
        mpq_t a2, sx, sy, cr, t1, t2, six;
        mpq_inits(a2, sx, sy, cr, t1, t2, six, (mpq_ptr)NULL);
        geom_signed_area2_q(&p, a2);
        if (mpq_sgn(a2) == 0) {
            mpq_clears(a2, sx, sy, cr, t1, t2, six, (mpq_ptr)NULL);
            geom_points_clear(&p);
            return NULL;  /* zero-area: documented deviation, decline */
        }
        for (size_t i = 0; i < p.n; i++) {
            size_t j = (i + 1) % p.n;
            mpq_mul(t1, p.qx[i], p.qy[j]);
            mpq_mul(t2, p.qx[j], p.qy[i]);
            mpq_sub(cr, t1, t2);                 /* cross_i */
            mpq_add(t1, p.qx[i], p.qx[j]);
            mpq_mul(t1, t1, cr);
            mpq_add(sx, sx, t1);
            mpq_add(t1, p.qy[i], p.qy[j]);
            mpq_mul(t1, t1, cr);
            mpq_add(sy, sy, t1);
        }
        mpq_set_ui(six, 3, 1);
        mpq_mul(six, six, a2);                   /* 6A == 3 * (2A) */
        mpq_div(sx, sx, six);
        mpq_div(sy, sy, six);
        out = geom_list2(geom_mpq_to_expr(sx), geom_mpq_to_expr(sy));
        mpq_clears(a2, sx, sy, cr, t1, t2, six, (mpq_ptr)NULL);
    } else {
        double a2 = geom_signed_area2_d(&p);
        if (a2 == 0.0) { geom_points_clear(&p); return NULL; }
        double sx = 0.0, sy = 0.0;
        for (size_t i = 0; i < p.n; i++) {
            size_t j = (i + 1) % p.n;
            double cr = p.dx[i] * p.dy[j] - p.dx[j] * p.dy[i];
            sx += (p.dx[i] + p.dx[j]) * cr;
            sy += (p.dy[i] + p.dy[j]) * cr;
        }
        out = geom_list2(expr_new_real(sx / (3.0 * a2)),
                         expr_new_real(sy / (3.0 * a2)));
    }
    geom_points_clear(&p);
    return out;
}

/* ------------------------------------------------------------------------- */
/* RegionMember                                                              */
/* ------------------------------------------------------------------------- */

/* Read RegionMember's second argument: a 2-element numeric List (or a rank-1,
 * length-2 NDArray). Fills both representations; *exact says whether qx/qy
 * are valid (caller must clear them only when 1). */
static int geom_read_query(const Expr* q, int* exact,
                           mpq_t qx, mpq_t qy, double* dx, double* dy) {
    if (is_ndarray(q)) {
        int n = 0;
        double* buf = NULL;
        if (!na_load_vector(q, false, &n, &buf)) return 0;
        if (n != 2) { free(buf); return 0; }
        *exact = 0;
        *dx = buf[0]; *dy = buf[1];
        free(buf);
        return 1;
    }
    if (!head_is(q, SYM_List) || q->data.function.arg_count != 2) return 0;
    const Expr* ex = q->data.function.args[0];
    const Expr* ey = q->data.function.args[1];
    if (geom_leaf_is_exact(ex) && geom_leaf_is_exact(ey)) {
        *exact = 1;
        mpq_init(qx); mpq_init(qy);
        geom_leaf_to_mpq(ex, qx);
        geom_leaf_to_mpq(ey, qy);
        *dx = mpq_get_d(qx); *dy = mpq_get_d(qy);
        return 1;
    }
    *exact = 0;
    return geom_leaf_to_double(ex, dx) && geom_leaf_to_double(ey, dy);
}

/* Boundary-inclusive point-in-polygon, exact path. */
static int geom_member_q(const GeomPoints* p, const mpq_t px, const mpq_t py) {
    int inside = 0;
    for (size_t i = 0; i < p->n; i++) {
        size_t j = (i + 1) % p->n;
        int cs = geom_cross_sign_q(p->qx[i], p->qy[i], p->qx[j], p->qy[j], px, py);
        if (cs == 0) {
            /* Collinear: on the segment iff within its bounding box. Written as
             * direct comparisons rather than `const mpq_t *lo = c ? &a : &b`:
             * mpq_t is an array type, so that conditional is a C99 constraint
             * violation (it only became legal in C23) and draws -Wpedantic
             * warnings that -Wall -Wextra do not show. */
            int in_x = (mpq_cmp(p->qx[i], px) <= 0 && mpq_cmp(px, p->qx[j]) <= 0) ||
                       (mpq_cmp(p->qx[j], px) <= 0 && mpq_cmp(px, p->qx[i]) <= 0);
            int in_y = (mpq_cmp(p->qy[i], py) <= 0 && mpq_cmp(py, p->qy[j]) <= 0) ||
                       (mpq_cmp(p->qy[j], py) <= 0 && mpq_cmp(py, p->qy[i]) <= 0);
            if (in_x && in_y) return 1;
        }
        int ai = mpq_cmp(p->qy[i], py) > 0;
        int bj = mpq_cmp(p->qy[j], py) > 0;
        if (ai != bj) {
            int s = mpq_cmp(p->qy[j], p->qy[i]);  /* sign of by - ay, nonzero */
            s = (s > 0) - (s < 0);
            if (cs == s) inside = !inside;
        }
    }
    return inside;
}

/* Boundary-inclusive point-in-polygon, machine path (exact IEEE compares). */
static int geom_member_d(const GeomPoints* p, double px, double py) {
    int inside = 0;
    for (size_t i = 0; i < p->n; i++) {
        size_t j = (i + 1) % p->n;
        int cs = geom_cross_sign_d(p->dx[i], p->dy[i], p->dx[j], p->dy[j], px, py);
        if (cs == 0) {
            double xmin = p->dx[i] < p->dx[j] ? p->dx[i] : p->dx[j];
            double xmax = p->dx[i] < p->dx[j] ? p->dx[j] : p->dx[i];
            double ymin = p->dy[i] < p->dy[j] ? p->dy[i] : p->dy[j];
            double ymax = p->dy[i] < p->dy[j] ? p->dy[j] : p->dy[i];
            if (xmin <= px && px <= xmax && ymin <= py && py <= ymax) return 1;
        }
        int ai = p->dy[i] > py;
        int bj = p->dy[j] > py;
        if (ai != bj) {
            double dyy = p->dy[j] - p->dy[i];
            int s = (dyy > 0.0) - (dyy < 0.0);
            if (cs == s) inside = !inside;
        }
    }
    return inside;
}

Expr* builtin_region_member(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    const Expr* poly = res->data.function.args[0];
    if (!head_is(poly, SYM_Polygon) || poly->data.function.arg_count != 1) return NULL;

    GeomPoints p;
    if (!geom_read_points(poly->data.function.args[0], &p)) return NULL;
    geom_dedup_consecutive(&p);
    if (p.n < 3) { geom_points_clear(&p); return NULL; }

    mpq_t qx, qy;
    double dx = 0.0, dy = 0.0;
    int q_exact = 0;
    if (!geom_read_query(res->data.function.args[1], &q_exact, qx, qy, &dx, &dy)) {
        geom_points_clear(&p);
        return NULL;
    }

    int member;
    if (p.exact && q_exact) {
        member = geom_member_q(&p, qx, qy);
    } else if (!p.d_finite) {
        /* Exact polygon whose coordinates do not fit a double, asked about with
         * a machine point: the machine predicate would compute NaN crosses and
         * answer confidently wrong. Decline instead. */
        if (q_exact) { mpq_clear(qx); mpq_clear(qy); }
        geom_points_clear(&p);
        return NULL;
    } else {
        member = geom_member_d(&p, dx, dy);
    }

    if (q_exact) { mpq_clear(qx); mpq_clear(qy); }
    geom_points_clear(&p);
    return expr_new_symbol(member ? SYM_True : SYM_False);
}

/* ------------------------------------------------------------------------- */
/* ConvexHullRegion                                                          */
/* ------------------------------------------------------------------------- */

/* Andrew's monotone chain over indices into a GeomPoints, exact and machine
 * variants. Both sort lexicographically (x, then y), dedup, then build the
 * lower and upper chains dropping collinear middle points (cross <= 0 pop),
 * yielding the hull CCW starting from the lexicographic minimum -- which is
 * exactly the vertex order WL's ConvexHullRegion prints (verified per-AC). */

static const GeomPoints* g_sort_pts;  /* qsort context (no qsort_r in C99) */

static int geom_idx_cmp_q(const void* a, const void* b) {
    size_t i = *(const size_t*)a, j = *(const size_t*)b;
    int c = mpq_cmp(g_sort_pts->qx[i], g_sort_pts->qx[j]);
    if (c) return c;
    return mpq_cmp(g_sort_pts->qy[i], g_sort_pts->qy[j]);
}

static int geom_idx_cmp_d(const void* a, const void* b) {
    size_t i = *(const size_t*)a, j = *(const size_t*)b;
    if (g_sort_pts->dx[i] < g_sort_pts->dx[j]) return -1;
    if (g_sort_pts->dx[i] > g_sort_pts->dx[j]) return 1;
    if (g_sort_pts->dy[i] < g_sort_pts->dy[j]) return -1;
    if (g_sort_pts->dy[i] > g_sort_pts->dy[j]) return 1;
    return 0;
}

static int geom_pt_equal(const GeomPoints* p, size_t i, size_t j) {
    if (p->exact)
        return mpq_equal(p->qx[i], p->qx[j]) && mpq_equal(p->qy[i], p->qy[j]);
    return p->dx[i] == p->dx[j] && p->dy[i] == p->dy[j];
}

static int geom_turn(const GeomPoints* p, size_t a, size_t b, size_t c) {
    if (p->exact)
        return geom_cross_sign_q(p->qx[a], p->qy[a], p->qx[b], p->qy[b],
                                 p->qx[c], p->qy[c]);
    return geom_cross_sign_d(p->dx[a], p->dy[a], p->dx[b], p->dy[b],
                             p->dx[c], p->dy[c]);
}

Expr* builtin_convex_hull_region(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;

    GeomPoints p;
    if (!geom_read_points(res->data.function.args[0], &p)) return NULL;

    size_t n = p.n;
    size_t* idx = (size_t*)malloc(n * sizeof(size_t));
    size_t* hull = (size_t*)malloc((2 * n + 1) * sizeof(size_t));
    if (!idx || !hull) { free(idx); free(hull); geom_points_clear(&p); return NULL; }
    for (size_t i = 0; i < n; i++) idx[i] = i;

    g_sort_pts = &p;
    qsort(idx, n, sizeof(size_t), p.exact ? geom_idx_cmp_q : geom_idx_cmp_d);
    g_sort_pts = NULL;

    /* Dedup (sorted order makes duplicates adjacent). */
    size_t m = 0;
    for (size_t i = 0; i < n; i++)
        if (m == 0 || !geom_pt_equal(&p, idx[m - 1], idx[i])) idx[m++] = idx[i];

    size_t k = 0;
    if (m == 1) {
        hull[k++] = idx[0];
    } else {
        for (size_t i = 0; i < m; i++) {            /* lower chain */
            while (k >= 2 && geom_turn(&p, hull[k - 2], hull[k - 1], idx[i]) <= 0) k--;
            hull[k++] = idx[i];
        }
        size_t lower = k + 1;
        for (size_t ii = m - 1; ii-- > 0; ) {       /* upper chain */
            while (k >= lower && geom_turn(&p, hull[k - 2], hull[k - 1], idx[ii]) <= 0) k--;
            hull[k++] = idx[ii];
        }
        k--;  /* last point repeats the first */
    }

    /* Build the result: Point / Line / Polygon by hull size. */
    Expr* out = NULL;
    Expr** verts = (Expr**)malloc((k ? k : 1) * sizeof(Expr*));
    if (!verts) { free(idx); free(hull); geom_points_clear(&p); return NULL; }
    for (size_t i = 0; i < k; i++) {
        size_t v = hull[i];
        Expr* cx = p.exact ? geom_mpq_to_expr(p.qx[v]) : expr_new_real(p.dx[v]);
        Expr* cy = p.exact ? geom_mpq_to_expr(p.qy[v]) : expr_new_real(p.dy[v]);
        verts[i] = geom_list2(cx, cy);
        if (!verts[i]) {  /* OOM: a NULL must never reach expr_new_function */
            for (size_t k = 0; k < i; k++) expr_free(verts[k]);
            free(verts); free(idx); free(hull);
            geom_points_clear(&p);
            return NULL;
        }
    }
    free(idx); free(hull);

    if (k == 1) {
        out = expr_new_function(expr_new_symbol(SYM_Point), verts, 1);
        free(verts);
    } else {
        Expr* vlist = expr_new_function(expr_new_symbol(SYM_List), verts, k);
        free(verts);
        Expr** warg = (Expr**)malloc(sizeof(Expr*));
        if (!warg) { expr_free(vlist); geom_points_clear(&p); return NULL; }
        warg[0] = vlist;
        out = expr_new_function(
            expr_new_symbol(k == 2 ? SYM_Line : SYM_Polygon), warg, 1);
        free(warg);
    }
    geom_points_clear(&p);
    return out;
}

/* ------------------------------------------------------------------------- */
/* Registration                                                              */
/* ------------------------------------------------------------------------- */

void geometry_init(void) {
    symtab_add_builtin("Area", builtin_area);
    symtab_add_builtin("Perimeter", builtin_perimeter);
    symtab_add_builtin("RegionCentroid", builtin_region_centroid);
    symtab_add_builtin("RegionMember", builtin_region_member);
    symtab_add_builtin("ConvexHullRegion", builtin_convex_hull_region);

    /* WL: Attributes[Area] == {Protected, ReadProtected}, same for all five.
     * Deliberately NOT Listable (these are structural, not element-wise). */
    static const char* heads[] = {
        "Area", "Perimeter", "RegionCentroid", "RegionMember", "ConvexHullRegion"
    };
    for (size_t i = 0; i < sizeof heads / sizeof heads[0]; i++)
        symtab_get_def(heads[i])->attributes |= (ATTR_PROTECTED | ATTR_READPROTECTED);

    symtab_set_docstring("Area",
        "Area[Polygon[{{x1, y1}, ...}]] gives the area of a simple 2D polygon. "
        "Exact (Integer/Rational) coordinates give an exact result; any Real "
        "coordinate gives a machine-precision result. A polygon with fewer "
        "than 3 distinct vertices has Undefined area.");
    symtab_set_docstring("Perimeter",
        "Perimeter[Polygon[{{x1, y1}, ...}]] gives the perimeter of a simple "
        "2D polygon: the sum of its edge lengths, including the closing edge. "
        "Exact coordinates give an exact (possibly symbolic, e.g. 2 + Sqrt[2]) "
        "result; Real coordinates give a machine-precision result.");
    symtab_set_docstring("RegionCentroid",
        "RegionCentroid[Polygon[{{x1, y1}, ...}]] gives the centroid {cx, cy} "
        "of a simple 2D polygon with nonzero area. Exact coordinates give an "
        "exact result; Real coordinates give a machine-precision result.");
    symtab_set_docstring("RegionMember",
        "RegionMember[Polygon[{{x1, y1}, ...}], {x, y}] gives True if the "
        "point lies inside or on the boundary of the simple 2D polygon, and "
        "False otherwise.");
    symtab_set_docstring("ConvexHullRegion",
        "ConvexHullRegion[{{x1, y1}, ...}] gives the convex hull of a set of "
        "2D points: a Polygon with the hull vertices in counterclockwise "
        "order, a Line for collinear input, or a Point for a single point.");
}
