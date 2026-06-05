/*
 * interp_mpfr.c
 *
 * Arbitrary-precision (MPFR) kernels for InterpolatingFunction / Interpolation.
 * A faithful mirror of the double kernels in interp.c: same algorithms
 * (windowed Newton, natural cubic spline, tensor-product Hermite with
 * finite-difference estimation of unsupplied mixed partials), carried out in
 * mpfr_t at a chosen working precision.  The whole translation unit compiles to
 * nothing when built without MPFR, like the linalg *_mpfr.c units.
 *
 * Entry point (declared in interp.c, used by interp_apply):
 *   Expr* interp_eval_mpfr(domain, table, m, call_args, ders, orders,
 *                          method, Ksupplied, bits);
 */

#ifdef USE_MPFR

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpfr.h>

#include "interp.h"
#include "expr.h"
#include "sym_names.h"
#include "numeric.h"

#define RND MPFR_RNDN

enum { METHOD_DEFAULT = 0, METHOD_SPLINE = 1, METHOD_HERMITE = 2 };

/* --- structural helpers (shared shape with interp.c) ------------------ */

static bool m_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_List;
}

/* Coerce a real numeric Expr into the (pre-initialised) mpfr_t out.  Rejects
 * non-real / non-numeric values. */
static bool nm_set(const Expr* e, mpfr_t out) {
    mpfr_t im; mpfr_init2(im, mpfr_get_prec(out));
    bool inexact;
    bool ok = get_approx_mpfr(e, out, im, &inexact);
    if (ok && !mpfr_zero_p(im)) ok = false;
    mpfr_clear(im);
    return ok;
}

static Expr* m_entry_value(Expr* entry) { return entry->data.function.args[1]; }

static Expr* m_entry_coord(Expr* entry, size_t m, size_t k) {
    Expr* point = entry->data.function.args[0];
    if (m_is_list(point) && point->data.function.arg_count == m)
        return point->data.function.args[k];
    if (m == 1) return point;
    return NULL;
}

/* --- parsed grid (mpfr abscissae) ------------------------------------- */

typedef struct {
    size_t   m;
    size_t*  nk;
    mpfr_t** grid;       /* grid[k][i] */
    size_t*  stride;
    size_t   total;
    Expr**   entryAt;
    long     prec;
    size_t   galloc;     /* mpfr_t's initialised per grid[k] (= npts) */
    bool*    periodic;   /* per-dimension periodicity */
    mpfr_t*  period;     /* per-dimension period (grid span); init'd, 0 if aperiodic */
} MFun;

static void mfun_free(MFun* f) {
    if (f->grid) {
        for (size_t k = 0; k < f->m; k++) {
            if (f->grid[k]) { for (size_t i = 0; i < f->galloc; i++) mpfr_clear(f->grid[k][i]); free(f->grid[k]); }
        }
        free(f->grid);
    }
    if (f->period) { for (size_t k = 0; k < f->m; k++) mpfr_clear(f->period[k]); free(f->period); }
    free(f->nk); free(f->stride); free(f->entryAt); free(f->periodic);
    memset(f, 0, sizeof(*f));
}

/* sorted-unique insert of value v (mpfr) into g[0..*len-1] */
static void grid_insert_m(mpfr_t* g, size_t* len, const mpfr_t v) {
    size_t i = 0;
    while (i < *len && mpfr_less_p(g[i], v)) i++;
    if (i < *len && mpfr_equal_p(g[i], v)) return;
    for (size_t j = *len; j > i; j--) mpfr_swap(g[j], g[j - 1]);
    mpfr_set(g[i], v, RND);
    (*len)++;
}

static size_t grid_index_m(mpfr_t* g, size_t len, const mpfr_t v) {
    size_t lo = 0, hi = len;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (mpfr_less_p(g[mid], v)) lo = mid + 1; else hi = mid;
    }
    return (lo < len && mpfr_equal_p(g[lo], v)) ? lo : len;
}

/* largest i with grid[i] <= x, clamped to [0, n-2] */
static size_t bracket_m(mpfr_t* xs, size_t n, const mpfr_t x) {
    if (mpfr_lessequal_p(x, xs[0])) return 0;
    if (mpfr_greaterequal_p(x, xs[n - 1])) return n - 2;
    size_t lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        size_t mid = lo + (hi - lo) / 2;
        if (mpfr_lessequal_p(xs[mid], x)) lo = mid; else hi = mid;
    }
    return lo;
}

static bool build_grid_m(Expr* domain, Expr* table, size_t m, const bool* periodic,
                         long prec, MFun* out) {
    (void)domain;
    memset(out, 0, sizeof(*out));
    out->prec = prec;
    if (!m_is_list(table)) return false;
    size_t npts = table->data.function.arg_count;
    if (npts < 2) return false;

    size_t*  nk     = calloc(m, sizeof(size_t));
    mpfr_t** grid   = calloc(m, sizeof(mpfr_t*));
    size_t*  stride = calloc(m, sizeof(size_t));
    Expr**   eat    = NULL;
    if (!nk || !grid || !stride) goto fail;

    mpfr_t v; mpfr_init2(v, prec);
    for (size_t k = 0; k < m; k++) {
        grid[k] = malloc(sizeof(mpfr_t) * npts);
        if (!grid[k]) { mpfr_clear(v); goto fail; }
        for (size_t i = 0; i < npts; i++) mpfr_init2(grid[k][i], prec);
        size_t len = 0;
        for (size_t i = 0; i < npts; i++) {
            Expr* e = table->data.function.args[i];
            if (!m_is_list(e) || e->data.function.arg_count < 2) { mpfr_clear(v); goto fail_grid; }
            Expr* ck = m_entry_coord(e, m, k);
            if (!ck || !nm_set(ck, v)) { mpfr_clear(v); goto fail_grid; }
            grid_insert_m(grid[k], &len, v);
        }
        nk[k] = len;
        if (len < 2) { mpfr_clear(v); goto fail_grid; }
    }

    size_t total = 1;
    for (size_t k = 0; k < m; k++) total *= nk[k];
    stride[m - 1] = 1;
    for (size_t k = m - 1; k-- > 0; ) stride[k] = stride[k + 1] * nk[k + 1];

    eat = calloc(total, sizeof(Expr*));
    if (!eat) { mpfr_clear(v); goto fail_grid; }
    for (size_t i = 0; i < npts; i++) {
        Expr* e = table->data.function.args[i];
        size_t flat = 0;
        for (size_t k = 0; k < m; k++) {
            nm_set(m_entry_coord(e, m, k), v);
            flat += grid_index_m(grid[k], nk[k], v) * stride[k];
        }
        eat[flat] = e;
    }
    mpfr_clear(v);
    for (size_t i = 0; i < total; i++) if (!eat[i]) goto fail_grid;

    bool*   per  = calloc(m, sizeof(bool));
    mpfr_t* perd = malloc(sizeof(mpfr_t) * m);
    for (size_t k = 0; k < m; k++) {
        mpfr_init2(perd[k], prec);
        per[k] = periodic && periodic[k];
        if (per[k]) mpfr_sub(perd[k], grid[k][nk[k] - 1], grid[k][0], RND);
        else mpfr_set_zero(perd[k], 1);
    }

    out->m = m; out->nk = nk; out->grid = grid; out->stride = stride;
    out->total = total; out->entryAt = eat; out->galloc = npts;
    out->periodic = per; out->period = perd;
    return true;

fail_grid:
    for (size_t k = 0; k < m; k++)
        if (grid[k]) { for (size_t i = 0; i < npts; i++) mpfr_clear(grid[k][i]); free(grid[k]); }
    free(grid); free(nk); free(stride); free(eat);
    return false;
fail:
    free(nk); free(grid); free(stride);
    return false;
}

/* --- 1-D kernels (mpfr) ----------------------------------------------- */

static void newton_eval_m(mpfr_t* xs, mpfr_t* ys, size_t w, const mpfr_t p,
                          size_t d, mpfr_t out, long prec) {
    mpfr_t* c = malloc(sizeof(mpfr_t) * w);
    for (size_t k = 0; k < w; k++) { mpfr_init2(c[k], prec); mpfr_set(c[k], ys[k], RND); }
    mpfr_t t1; mpfr_init2(t1, prec);
    for (size_t j = 1; j < w; j++)
        for (size_t k = w - 1; k >= j; k--) {
            mpfr_sub(c[k], c[k], c[k - 1], RND);
            mpfr_sub(t1, xs[k], xs[k - j], RND);
            mpfr_div(c[k], c[k], t1, RND);
            if (k == j) break;
        }
    mpfr_t* val = malloc(sizeof(mpfr_t) * (d + 1));
    for (size_t i = 0; i <= d; i++) { mpfr_init2(val[i], prec); mpfr_set_zero(val[i], 1); }
    mpfr_set(val[0], c[w - 1], RND);
    mpfr_t dx; mpfr_init2(dx, prec);
    for (size_t ii = w - 1; ii-- > 0; ) {
        mpfr_sub(dx, p, xs[ii], RND);
        for (size_t t = d; t >= 1; t--) {
            mpfr_mul(val[t], val[t], dx, RND);
            mpfr_add(val[t], val[t], val[t - 1], RND);
            if (t == 1) break;
        }
        mpfr_mul(val[0], val[0], dx, RND);
        mpfr_add(val[0], val[0], c[ii], RND);
    }
    unsigned long f = 1;
    for (size_t t = 2; t <= d; t++) f *= (unsigned long)t;
    mpfr_mul_ui(out, val[d], f, RND);
    for (size_t k = 0; k < w; k++) mpfr_clear(c[k]);
    for (size_t i = 0; i <= d; i++) mpfr_clear(val[i]);
    mpfr_clear(t1); mpfr_clear(dx); free(c); free(val);
}

static void spline_eval_m(mpfr_t* xs, mpfr_t* ys, size_t n, const mpfr_t p,
                          size_t d, mpfr_t out, long prec) {
    mpfr_t* M = malloc(sizeof(mpfr_t) * n);
    for (size_t i = 0; i < n; i++) { mpfr_init2(M[i], prec); mpfr_set_zero(M[i], 1); }
    if (n >= 3) {
        size_t ni = n - 2;
        mpfr_t* cc = malloc(sizeof(mpfr_t) * ni);
        mpfr_t* rr = malloc(sizeof(mpfr_t) * ni);
        for (size_t j = 0; j < ni; j++) { mpfr_init2(cc[j], prec); mpfr_init2(rr[j], prec); }
        mpfr_t h0, h1, diag, sub, sup, d0, d1, rhs, denom, tmp;
        mpfr_inits2(prec, h0, h1, diag, sub, sup, d0, d1, rhs, denom, tmp, (mpfr_ptr)0);
        for (size_t j = 0; j < ni; j++) {
            size_t i = j + 1;
            mpfr_sub(h0, xs[i], xs[i - 1], RND);
            mpfr_sub(h1, xs[i + 1], xs[i], RND);
            mpfr_add(diag, h0, h1, RND); mpfr_mul_ui(diag, diag, 2, RND);
            mpfr_set(sub, h0, RND); mpfr_set(sup, h1, RND);
            mpfr_sub(d0, ys[i], ys[i - 1], RND); mpfr_div(d0, d0, h0, RND);
            mpfr_sub(d1, ys[i + 1], ys[i], RND); mpfr_div(d1, d1, h1, RND);
            mpfr_sub(rhs, d1, d0, RND); mpfr_mul_ui(rhs, rhs, 6, RND);
            if (j == 0) { mpfr_div(cc[j], sup, diag, RND); mpfr_div(rr[j], rhs, diag, RND); }
            else {
                mpfr_mul(tmp, sub, cc[j - 1], RND); mpfr_sub(denom, diag, tmp, RND);
                mpfr_div(cc[j], sup, denom, RND);
                mpfr_mul(tmp, sub, rr[j - 1], RND); mpfr_sub(tmp, rhs, tmp, RND);
                mpfr_div(rr[j], tmp, denom, RND);
            }
        }
        mpfr_set(M[n - 2], rr[ni - 1], RND);
        for (size_t j = ni - 1; j-- > 0; ) {
            mpfr_mul(tmp, cc[j], M[j + 2], RND);
            mpfr_sub(M[j + 1], rr[j], tmp, RND);
        }
        for (size_t j = 0; j < ni; j++) { mpfr_clear(cc[j]); mpfr_clear(rr[j]); }
        free(cc); free(rr);
        mpfr_clears(h0, h1, diag, sub, sup, d0, d1, rhs, denom, tmp, (mpfr_ptr)0);
    }
    size_t i = bracket_m(xs, n, p);
    mpfr_t h, A, B, Mi, Mi1, ci, ci1, t1, t2;
    mpfr_inits2(prec, h, A, B, Mi, Mi1, ci, ci1, t1, t2, (mpfr_ptr)0);
    mpfr_sub(h, xs[i + 1], xs[i], RND);
    mpfr_sub(A, xs[i + 1], p, RND);
    mpfr_sub(B, p, xs[i], RND);
    mpfr_set(Mi, M[i], RND); mpfr_set(Mi1, M[i + 1], RND);
    /* ci  = ys[i]/h   - Mi *h/6 ; ci1 = ys[i+1]/h - Mi1*h/6 */
    mpfr_div(ci, ys[i], h, RND);     mpfr_mul(t1, Mi, h, RND);  mpfr_div_ui(t1, t1, 6, RND); mpfr_sub(ci, ci, t1, RND);
    mpfr_div(ci1, ys[i + 1], h, RND); mpfr_mul(t1, Mi1, h, RND); mpfr_div_ui(t1, t1, 6, RND); mpfr_sub(ci1, ci1, t1, RND);
    if (d == 0) {
        /* Mi*A^3/(6h) + Mi1*B^3/(6h) + ci*A + ci1*B */
        mpfr_pow_ui(t1, A, 3, RND); mpfr_mul(t1, t1, Mi, RND); mpfr_div(t1, t1, h, RND); mpfr_div_ui(t1, t1, 6, RND); mpfr_set(out, t1, RND);
        mpfr_pow_ui(t1, B, 3, RND); mpfr_mul(t1, t1, Mi1, RND); mpfr_div(t1, t1, h, RND); mpfr_div_ui(t1, t1, 6, RND); mpfr_add(out, out, t1, RND);
        mpfr_mul(t1, ci, A, RND); mpfr_add(out, out, t1, RND);
        mpfr_mul(t1, ci1, B, RND); mpfr_add(out, out, t1, RND);
    } else if (d == 1) {
        /* -Mi*A^2/(2h) + Mi1*B^2/(2h) - ci + ci1 */
        mpfr_mul(t1, A, A, RND); mpfr_mul(t1, t1, Mi, RND); mpfr_div(t1, t1, h, RND); mpfr_div_ui(t1, t1, 2, RND); mpfr_neg(out, t1, RND);
        mpfr_mul(t1, B, B, RND); mpfr_mul(t1, t1, Mi1, RND); mpfr_div(t1, t1, h, RND); mpfr_div_ui(t1, t1, 2, RND); mpfr_add(out, out, t1, RND);
        mpfr_sub(out, out, ci, RND); mpfr_add(out, out, ci1, RND);
    } else if (d == 2) {
        /* Mi*A/h + Mi1*B/h */
        mpfr_mul(t1, Mi, A, RND); mpfr_div(t1, t1, h, RND); mpfr_set(out, t1, RND);
        mpfr_mul(t1, Mi1, B, RND); mpfr_div(t1, t1, h, RND); mpfr_add(out, out, t1, RND);
    } else if (d == 3) {
        mpfr_sub(t1, Mi1, Mi, RND); mpfr_div(out, t1, h, RND);
    } else {
        mpfr_set_zero(out, 1);
    }
    for (size_t k = 0; k < n; k++) mpfr_clear(M[k]);
    free(M);
    mpfr_clears(h, A, B, Mi, Mi1, ci, ci1, t1, t2, (mpfr_ptr)0);
}

/* Periodic (cyclic) cubic spline; mirror of spline_eval_periodic.  xs[0..n-1]
 * distinct nodes, period P (x_n = x_0 + P), p reduced into [x_0, x_0+P). */
static void spline_eval_periodic_m(mpfr_t* xs, size_t n, mpfr_t* ys, const mpfr_t P,
                                   const mpfr_t p, size_t d, mpfr_t out, long prec) {
    mpfr_t* h = malloc(sizeof(mpfr_t) * n);
    for (size_t i = 0; i < n; i++) mpfr_init2(h[i], prec);
    for (size_t i = 0; i + 1 < n; i++) mpfr_sub(h[i], xs[i + 1], xs[i], RND);
    mpfr_add(h[n - 1], xs[0], P, RND); mpfr_sub(h[n - 1], h[n - 1], xs[n - 1], RND);

    mpfr_t* M = malloc(sizeof(mpfr_t) * n);
    for (size_t i = 0; i < n; i++) { mpfr_init2(M[i], prec); mpfr_set_zero(M[i], 1); }

    if (n >= 3) {
        mpfr_t* a = malloc(sizeof(mpfr_t) * n);
        mpfr_t* b = malloc(sizeof(mpfr_t) * n);
        mpfr_t* c = malloc(sizeof(mpfr_t) * n);
        mpfr_t* r = malloc(sizeof(mpfr_t) * n);
        mpfr_t* bb = malloc(sizeof(mpfr_t) * n);
        mpfr_t* u = malloc(sizeof(mpfr_t) * n);
        mpfr_t* cp = malloc(sizeof(mpfr_t) * n);
        mpfr_t* y1 = malloc(sizeof(mpfr_t) * n);
        mpfr_t* z1 = malloc(sizeof(mpfr_t) * n);
        for (size_t i = 0; i < n; i++)
            mpfr_inits2(prec, a[i], b[i], c[i], r[i], bb[i], u[i], cp[i], y1[i], z1[i], (mpfr_ptr)0);
        mpfr_t dn, dp, tmp; mpfr_inits2(prec, dn, dp, tmp, (mpfr_ptr)0);
        for (size_t i = 0; i < n; i++) {
            mpfr_t* hm = &h[(i + n - 1) % n]; mpfr_t* hi = &h[i];
            mpfr_set(a[i], *hm, RND); mpfr_set(c[i], *hi, RND);
            mpfr_add(b[i], *hm, *hi, RND); mpfr_mul_ui(b[i], b[i], 2, RND);
            mpfr_sub(dn, ys[(i + 1) % n], ys[i], RND); mpfr_div(dn, dn, *hi, RND);
            mpfr_sub(dp, ys[i], ys[(i + n - 1) % n], RND); mpfr_div(dp, dp, *hm, RND);
            mpfr_sub(r[i], dn, dp, RND); mpfr_mul_ui(r[i], r[i], 6, RND);
        }
        /* Sherman-Morrison setup: alpha=a[0], beta=c[n-1], gamma=-b[0]. */
        mpfr_t alpha, beta, gamma, ab;
        mpfr_inits2(prec, alpha, beta, gamma, ab, (mpfr_ptr)0);
        mpfr_set(alpha, a[0], RND); mpfr_set(beta, c[n - 1], RND); mpfr_neg(gamma, b[0], RND);
        for (size_t i = 0; i < n; i++) { mpfr_set(bb[i], b[i], RND); mpfr_set_zero(u[i], 1); }
        mpfr_sub(bb[0], b[0], gamma, RND);
        mpfr_mul(ab, alpha, beta, RND); mpfr_div(ab, ab, gamma, RND);
        mpfr_sub(bb[n - 1], b[n - 1], ab, RND);
        mpfr_set(u[0], gamma, RND); mpfr_set(u[n - 1], beta, RND);
        /* Thomas solve T y1 = r and T z1 = u. */
        mpfr_div(cp[0], c[0], bb[0], RND); mpfr_div(y1[0], r[0], bb[0], RND); mpfr_div(z1[0], u[0], bb[0], RND);
        for (size_t i = 1; i < n; i++) {
            mpfr_mul(tmp, a[i], cp[i - 1], RND); mpfr_sub(tmp, bb[i], tmp, RND);  /* mden */
            mpfr_div(cp[i], c[i], tmp, RND);
            mpfr_t t2b; mpfr_init2(t2b, prec);
            mpfr_mul(t2b, a[i], y1[i - 1], RND); mpfr_sub(t2b, r[i], t2b, RND); mpfr_div(y1[i], t2b, tmp, RND);
            mpfr_mul(t2b, a[i], z1[i - 1], RND); mpfr_sub(t2b, u[i], t2b, RND); mpfr_div(z1[i], t2b, tmp, RND);
            mpfr_clear(t2b);
        }
        for (size_t i = n - 1; i-- > 0; ) {
            mpfr_mul(tmp, cp[i], y1[i + 1], RND); mpfr_sub(y1[i], y1[i], tmp, RND);
            mpfr_mul(tmp, cp[i], z1[i + 1], RND); mpfr_sub(z1[i], z1[i], tmp, RND);
        }
        /* x = y1 - (v.y1)/(1+v.z1) z1, v = e_0 + (alpha/gamma) e_{n-1}. */
        mpfr_t ag, vy, vz, fac; mpfr_inits2(prec, ag, vy, vz, fac, (mpfr_ptr)0);
        mpfr_div(ag, alpha, gamma, RND);
        mpfr_mul(vy, ag, y1[n - 1], RND); mpfr_add(vy, vy, y1[0], RND);
        mpfr_mul(vz, ag, z1[n - 1], RND); mpfr_add(vz, vz, z1[0], RND);
        mpfr_set_ui(fac, 1, RND); mpfr_add(fac, fac, vz, RND);   /* 1 + vz */
        mpfr_div(fac, vy, fac, RND);
        for (size_t i = 0; i < n; i++) { mpfr_mul(tmp, fac, z1[i], RND); mpfr_sub(M[i], y1[i], tmp, RND); }
        for (size_t i = 0; i < n; i++)
            mpfr_clears(a[i], b[i], c[i], r[i], bb[i], u[i], cp[i], y1[i], z1[i], (mpfr_ptr)0);
        free(a); free(b); free(c); free(r); free(bb); free(u); free(cp); free(y1); free(z1);
        mpfr_clears(dn, dp, tmp, alpha, beta, gamma, ab, ag, vy, vz, fac, (mpfr_ptr)0);
    }

    /* Locate interval i (0..n-1); interval n-1 is the wrap. */
    size_t i = 0;
    if (mpfr_greaterequal_p(p, xs[n - 1])) i = n - 1;
    else { while (i + 1 < n && mpfr_lessequal_p(xs[i + 1], p)) i++; }
    mpfr_t xR, yR, MR, hh, A, B, cL, cR, t1, t2;
    mpfr_inits2(prec, xR, yR, MR, hh, A, B, cL, cR, t1, t2, (mpfr_ptr)0);
    if (i + 1 < n) { mpfr_set(xR, xs[i + 1], RND); mpfr_set(yR, ys[i + 1], RND); mpfr_set(MR, M[i + 1], RND); }
    else { mpfr_add(xR, xs[0], P, RND); mpfr_set(yR, ys[0], RND); mpfr_set(MR, M[0], RND); }
    mpfr_sub(hh, xR, xs[i], RND);
    mpfr_sub(A, xR, p, RND); mpfr_sub(B, p, xs[i], RND);
    mpfr_div(cL, ys[i], hh, RND); mpfr_mul(t1, M[i], hh, RND); mpfr_div_ui(t1, t1, 6, RND); mpfr_sub(cL, cL, t1, RND);
    mpfr_div(cR, yR, hh, RND);    mpfr_mul(t1, MR, hh, RND);   mpfr_div_ui(t1, t1, 6, RND); mpfr_sub(cR, cR, t1, RND);
    if (d == 0) {
        mpfr_pow_ui(t1, A, 3, RND); mpfr_mul(t1, t1, M[i], RND); mpfr_div(t1, t1, hh, RND); mpfr_div_ui(t1, t1, 6, RND); mpfr_set(out, t1, RND);
        mpfr_pow_ui(t1, B, 3, RND); mpfr_mul(t1, t1, MR, RND); mpfr_div(t1, t1, hh, RND); mpfr_div_ui(t1, t1, 6, RND); mpfr_add(out, out, t1, RND);
        mpfr_mul(t1, cL, A, RND); mpfr_add(out, out, t1, RND);
        mpfr_mul(t1, cR, B, RND); mpfr_add(out, out, t1, RND);
    } else if (d == 1) {
        mpfr_mul(t1, A, A, RND); mpfr_mul(t1, t1, M[i], RND); mpfr_div(t1, t1, hh, RND); mpfr_div_ui(t1, t1, 2, RND); mpfr_neg(out, t1, RND);
        mpfr_mul(t1, B, B, RND); mpfr_mul(t1, t1, MR, RND); mpfr_div(t1, t1, hh, RND); mpfr_div_ui(t1, t1, 2, RND); mpfr_add(out, out, t1, RND);
        mpfr_sub(out, out, cL, RND); mpfr_add(out, out, cR, RND);
    } else if (d == 2) {
        mpfr_mul(t1, M[i], A, RND); mpfr_div(t1, t1, hh, RND); mpfr_set(out, t1, RND);
        mpfr_mul(t1, MR, B, RND); mpfr_div(t1, t1, hh, RND); mpfr_add(out, out, t1, RND);
    } else if (d == 3) {
        mpfr_sub(t1, MR, M[i], RND); mpfr_div(out, t1, hh, RND);
    } else mpfr_set_zero(out, 1);

    for (size_t k = 0; k < n; k++) { mpfr_clear(h[k]); mpfr_clear(M[k]); }
    free(h); free(M);
    mpfr_clears(xR, yR, MR, hh, A, B, cL, cR, t1, t2, (mpfr_ptr)0);
}

/* --- Hermite tensor machinery (mpfr) ---------------------------------- */

static size_t ipow_sz(size_t base, size_t e) { size_t r = 1; while (e--) r *= base; return r; }
static void decode_mi(size_t idx, size_t m, int k, int* a) {
    for (size_t d = m; d-- > 0; ) { a[d] = (int)(idx % (size_t)(k + 1)); idx /= (size_t)(k + 1); }
}
static int mi_sum(const int* a, size_t m) { int s = 0; for (size_t d = 0; d < m; d++) s += a[d]; return s; }

/* falling factorial p*(p-1)*...*(p-g+1) as an mpfr at out's precision */
static void falling_m(int p, int g, mpfr_t out) {
    mpfr_set_ui(out, 1, RND);
    for (int i = 0; i < g; i++) mpfr_mul_si(out, out, p - i, RND);
}

/*
 * build_basis_m:
 *   The 2(k+1) one-dimensional Hermite basis polynomials on [0,1] at `prec`.
 *   basis[(e*(k+1)+j)*(2k+2) + p] = t^p coefficient of phi_{e,j}.  Mirrors
 *   build_basis() in interp.c.  Caller mpfr_clear's each coeff and frees.
 */
static mpfr_t* build_basis_m(int k, long prec) {
    int D = 2 * k + 2, K1 = k + 1;
    size_t ncoef = (size_t)2 * K1 * D;
    mpfr_t* basis = malloc(sizeof(mpfr_t) * ncoef);
    for (size_t i = 0; i < ncoef; i++) { mpfr_init2(basis[i], prec); mpfr_set_zero(basis[i], 1); }
    mpfr_t fac, fl, t, g1;
    mpfr_inits2(prec, fac, fl, t, g1, (mpfr_ptr)0);
    /* A[K1][K1], rhs[K1] as mpfr */
    mpfr_t A[16][16], rhs[16];
    for (int i = 0; i < K1; i++) { mpfr_init2(rhs[i], prec); for (int q = 0; q < K1; q++) mpfr_init2(A[i][q], prec); }

    for (int e = 0; e < 2; e++) {
        for (int j = 0; j < K1; j++) {
            mpfr_t* c = basis + ((size_t)(e * K1 + j)) * D;
            for (int i = 0; i < K1; i++) {
                if (e == 0 && i == j) { mpfr_set_ui(fac, 1, RND); for (int q = 2; q <= i; q++) mpfr_mul_ui(fac, fac, q, RND);
                                        mpfr_ui_div(c[i], 1, fac, RND); }
                else mpfr_set_zero(c[i], 1);
            }
            for (int i = 0; i < K1; i++) {
                mpfr_set_zero(t, 1);
                for (int p = i; p <= k; p++) { falling_m(p, i, fl); mpfr_mul(fl, fl, c[p], RND); mpfr_add(t, t, fl, RND); }
                mpfr_set_si(g1, (e == 1 && i == j) ? 1 : 0, RND);
                mpfr_sub(rhs[i], g1, t, RND);
                for (int q = 0; q < K1; q++) { int p = k + 1 + q; if (p >= i) falling_m(p, i, A[i][q]); else mpfr_set_zero(A[i][q], 1); }
            }
            /* Gaussian elimination with partial pivoting */
            for (int col = 0; col < K1; col++) {
                int piv = col;
                for (int r = col + 1; r < K1; r++) if (mpfr_cmpabs(A[r][col], A[piv][col]) > 0) piv = r;
                if (piv != col) {
                    for (int cc = 0; cc < K1; cc++) mpfr_swap(A[col][cc], A[piv][cc]);
                    mpfr_swap(rhs[col], rhs[piv]);
                }
                for (int r = 0; r < K1; r++) if (r != col) {
                    mpfr_div(fac, A[r][col], A[col][col], RND);
                    for (int cc = col; cc < K1; cc++) { mpfr_mul(fl, fac, A[col][cc], RND); mpfr_sub(A[r][cc], A[r][cc], fl, RND); }
                    mpfr_mul(fl, fac, rhs[col], RND); mpfr_sub(rhs[r], rhs[r], fl, RND);
                }
            }
            for (int q = 0; q < K1; q++) mpfr_div(c[k + 1 + q], rhs[q], A[q][q], RND);
        }
    }
    for (int i = 0; i < K1; i++) { mpfr_clear(rhs[i]); for (int q = 0; q < K1; q++) mpfr_clear(A[i][q]); }
    mpfr_clears(fac, fl, t, g1, (mpfr_ptr)0);
    return basis;
}

/* g-th derivative at t of polynomial coeffs c[0..deg], into out */
static void phi_eval_m(mpfr_t* c, int deg, int g, const mpfr_t t, mpfr_t out, long prec) {
    mpfr_t fl; mpfr_init2(fl, prec);
    mpfr_set_zero(out, 1);
    for (int p = deg; p >= g; p--) {
        mpfr_mul(out, out, t, RND);
        falling_m(p, g, fl); mpfr_mul(fl, fl, c[p], RND);
        mpfr_add(out, out, fl, RND);
    }
    mpfr_clear(fl);
}

static Expr* descend_cpath_m(Expr* t, const int* cpath, int vrank) {
    for (int d = 0; d < vrank; d++) {
        if (!m_is_list(t) || (int)t->data.function.arg_count <= cpath[d]) return NULL;
        t = t->data.function.args[cpath[d]];
    }
    return t;
}

static bool value_component_m(Expr* entry, const int* cpath, int vrank, mpfr_t out) {
    Expr* t = descend_cpath_m(entry->data.function.args[1], cpath, vrank);
    return t && nm_set(t, out);
}

/* Value shape from the first node value; rank 0 for a scalar. */
static void value_shape_m(Expr* e, size_t* vshape, int* vrank) {
    int r = 0;
    while (m_is_list(e) && r < 16) { vshape[r++] = e->data.function.arg_count; e = e->data.function.args[0]; }
    *vrank = r;
}

/* D^a f_component at node `entry` (a has total order s); from the supplied
 * tensor, after selecting the value-array component `cpath`. */
static bool extract_supplied_m(Expr* entry, size_t m, const int* cpath, int vrank,
                               const int* a, int s, mpfr_t out) {
    if ((size_t)(1 + s) >= entry->data.function.arg_count) return false;
    Expr* t = descend_cpath_m(entry->data.function.args[1 + s], cpath, vrank);
    if (!t) return false;
    if (m == 1) {
        while (m_is_list(t) && t->data.function.arg_count == 1) t = t->data.function.args[0];
        return nm_set(t, out);
    }
    for (size_t d = 0; d < m; d++)
        for (int r = 0; r < a[d]; r++) {
            if (!m_is_list(t) || t->data.function.arg_count != m) return false;
            t = t->data.function.args[d];
        }
    return nm_set(t, out);
}

/* estimate d/dx_d of field F over the grid, into G (3-point / one-sided;
 * cyclic central differences on a periodic dimension) */
static void fd_apply_m(const MFun* g, size_t d, mpfr_t* F, mpfr_t* G, long prec) {
    size_t nd = g->nk[d], st = g->stride[d];
    mpfr_t* xs = g->grid[d];
    mpfr_t h0, h1, d0, d1, t1, t2;
    mpfr_inits2(prec, h0, h1, d0, d1, t1, t2, (mpfr_ptr)0);
    if (g->periodic && g->periodic[d] && nd >= 3) {
        size_t ndist = nd - 1;
        for (size_t n = 0; n < g->total; n++) {
            size_t id = (n / st) % nd;
            size_t base = n - id * st;
            size_t idm = (id == nd - 1) ? 0 : id;
            size_t L = (idm + ndist - 1) % ndist, R = (idm + 1) % ndist;
            mpfr_sub(h0, xs[idm], xs[L], RND); if (idm == 0) mpfr_add(h0, h0, g->period[d], RND);
            mpfr_sub(h1, xs[R], xs[idm], RND); if (idm == ndist - 1) mpfr_add(h1, h1, g->period[d], RND);
            mpfr_sub(d0, F[base + idm * st], F[base + L * st], RND); mpfr_div(d0, d0, h0, RND);
            mpfr_sub(d1, F[base + R * st], F[base + idm * st], RND); mpfr_div(d1, d1, h1, RND);
            mpfr_mul(t1, h1, d0, RND); mpfr_mul(t2, h0, d1, RND); mpfr_add(t1, t1, t2, RND);
            mpfr_add(t2, h0, h1, RND); mpfr_div(G[n], t1, t2, RND);
        }
        mpfr_clears(h0, h1, d0, d1, t1, t2, (mpfr_ptr)0);
        return;
    }
    for (size_t n = 0; n < g->total; n++) {
        size_t id = (n / st) % nd;
        size_t base = n - id * st;
        if (nd == 1) { mpfr_set_zero(G[n], 1); continue; }
        if (nd == 2) {
            mpfr_sub(t1, F[base + st], F[base], RND);
            mpfr_sub(t2, xs[1], xs[0], RND);
            mpfr_div(G[n], t1, t2, RND);
            continue;
        }
        if (id == 0) {
            mpfr_sub(h0, xs[1], xs[0], RND); mpfr_sub(h1, xs[2], xs[1], RND);
            mpfr_sub(d0, F[base + st], F[base], RND); mpfr_div(d0, d0, h0, RND);
            mpfr_sub(d1, F[base + 2 * st], F[base + st], RND); mpfr_div(d1, d1, h1, RND);
            /* ((2h0+h1)d0 - h0 d1)/(h0+h1) */
            mpfr_mul_ui(t1, h0, 2, RND); mpfr_add(t1, t1, h1, RND); mpfr_mul(t1, t1, d0, RND);
            mpfr_mul(t2, h0, d1, RND); mpfr_sub(t1, t1, t2, RND);
            mpfr_add(t2, h0, h1, RND); mpfr_div(G[n], t1, t2, RND);
        } else if (id == nd - 1) {
            mpfr_sub(h0, xs[nd - 2], xs[nd - 3], RND); mpfr_sub(h1, xs[nd - 1], xs[nd - 2], RND);
            mpfr_sub(d0, F[base + (nd - 2) * st], F[base + (nd - 3) * st], RND); mpfr_div(d0, d0, h0, RND);
            mpfr_sub(d1, F[base + (nd - 1) * st], F[base + (nd - 2) * st], RND); mpfr_div(d1, d1, h1, RND);
            /* ((2h1+h0)d1 - h1 d0)/(h0+h1) */
            mpfr_mul_ui(t1, h1, 2, RND); mpfr_add(t1, t1, h0, RND); mpfr_mul(t1, t1, d1, RND);
            mpfr_mul(t2, h1, d0, RND); mpfr_sub(t1, t1, t2, RND);
            mpfr_add(t2, h0, h1, RND); mpfr_div(G[n], t1, t2, RND);
        } else {
            mpfr_sub(h0, xs[id], xs[id - 1], RND); mpfr_sub(h1, xs[id + 1], xs[id], RND);
            mpfr_sub(d0, F[base + id * st], F[base + (id - 1) * st], RND); mpfr_div(d0, d0, h0, RND);
            mpfr_sub(d1, F[base + (id + 1) * st], F[base + id * st], RND); mpfr_div(d1, d1, h1, RND);
            /* (h1 d0 + h0 d1)/(h0+h1) */
            mpfr_mul(t1, h1, d0, RND); mpfr_mul(t2, h0, d1, RND); mpfr_add(t1, t1, t2, RND);
            mpfr_add(t2, h0, h1, RND); mpfr_div(G[n], t1, t2, RND);
        }
    }
    mpfr_clears(h0, h1, d0, d1, t1, t2, (mpfr_ptr)0);
}

/* full mixed-derivative tensor T[node*dcount+idx], idx over {0..k}^m */
static mpfr_t* build_T_m(const MFun* g, size_t m, int K, int k, size_t dcount,
                         const int* cpath, int vrank, long prec) {
    mpfr_t* T = malloc(sizeof(mpfr_t) * g->total * dcount);
    for (size_t i = 0; i < g->total * dcount; i++) { mpfr_init2(T[i], prec); mpfr_set_zero(T[i], 1); }
    mpfr_t* F = malloc(sizeof(mpfr_t) * g->total);
    mpfr_t* Gt = malloc(sizeof(mpfr_t) * g->total);
    for (size_t i = 0; i < g->total; i++) { mpfr_init2(F[i], prec); mpfr_init2(Gt[i], prec); }
    int* a = malloc(sizeof(int) * m);
    int* b = malloc(sizeof(int) * m);
    bool ok = true;
    for (size_t idx = 0; ok && idx < dcount; idx++) {
        decode_mi(idx, m, k, a);
        int s = mi_sum(a, m);
        if (s <= K) {
            for (size_t n = 0; n < g->total; n++)
                if (!extract_supplied_m(g->entryAt[n], m, cpath, vrank, a, s, T[n * dcount + idx])) { ok = false; break; }
        } else {
            int rem = K;
            for (size_t d = 0; d < m; d++) { int bd = a[d] < rem ? a[d] : rem; b[d] = bd; rem -= bd; }
            for (size_t n = 0; n < g->total; n++)
                if (!extract_supplied_m(g->entryAt[n], m, cpath, vrank, b, K, F[n])) { ok = false; break; }
            if (!ok) break;
            for (size_t d = 0; d < m; d++) {
                int e = a[d] - b[d];
                for (int r = 0; r < e; r++) { fd_apply_m(g, d, F, Gt, prec); mpfr_t* sw = F; F = Gt; Gt = sw; }
            }
            for (size_t n = 0; n < g->total; n++) mpfr_set(T[n * dcount + idx], F[n], RND);
        }
    }
    for (size_t i = 0; i < g->total; i++) { mpfr_clear(F[i]); mpfr_clear(Gt[i]); }
    free(F); free(Gt); free(a); free(b);
    if (!ok) { for (size_t i = 0; i < g->total * dcount; i++) mpfr_clear(T[i]); free(T); return NULL; }
    return T;
}

static void hermite_tensor_eval_m(const MFun* g, mpfr_t* T, size_t dcount, size_t m,
                                  int k, mpfr_t* basis, int D, mpfr_t* p,
                                  const int* ders, mpfr_t out, long prec) {
    int* lo = malloc(sizeof(int) * m);
    int* a = malloc(sizeof(int) * m);
    int* c = malloc(sizeof(int) * m);
    mpfr_t* h = malloc(sizeof(mpfr_t) * m);
    mpfr_t* t = malloc(sizeof(mpfr_t) * m);
    size_t fcount = m * 2 * (size_t)(k + 1);
    mpfr_t* factor = malloc(sizeof(mpfr_t) * fcount);
    for (size_t i = 0; i < m; i++) { mpfr_init2(h[i], prec); mpfr_init2(t[i], prec); }
    for (size_t i = 0; i < fcount; i++) mpfr_init2(factor[i], prec);
    mpfr_t phi, hp, term, tmp;
    mpfr_inits2(prec, phi, hp, term, tmp, (mpfr_ptr)0);

    for (size_t d = 0; d < m; d++) {
        size_t i = bracket_m(g->grid[d], g->nk[d], p[d]);
        lo[d] = (int)i;
        mpfr_sub(h[d], g->grid[d][i + 1], g->grid[d][i], RND);
        mpfr_sub(t[d], p[d], g->grid[d][i], RND); mpfr_div(t[d], t[d], h[d], RND);
        for (int cd = 0; cd < 2; cd++)
            for (int ad = 0; ad <= k; ad++) {
                mpfr_t* coeffs = basis + ((size_t)(cd * (k + 1) + ad)) * D;
                phi_eval_m(coeffs, 2 * k + 1, ders[d], t[d], phi, prec);
                /* hp = h[d]^(ad - ders[d]) */
                int ex = ad - ders[d];
                if (ex >= 0) mpfr_pow_ui(hp, h[d], (unsigned long)ex, RND);
                else { mpfr_pow_ui(hp, h[d], (unsigned long)(-ex), RND); mpfr_ui_div(hp, 1, hp, RND); }
                mpfr_mul(factor[(d * 2 + cd) * (size_t)(k + 1) + ad], phi, hp, RND);
            }
    }

    mpfr_set_zero(out, 1);
    size_t ncorner = (size_t)1 << m;
    for (size_t cc = 0; cc < ncorner; cc++) {
        for (size_t d = 0; d < m; d++) c[d] = (int)((cc >> d) & 1);
        size_t node = 0;
        for (size_t d = 0; d < m; d++) node += (size_t)(lo[d] + c[d]) * g->stride[d];
        for (size_t idx = 0; idx < dcount; idx++) {
            if (mpfr_zero_p(T[node * dcount + idx])) continue;
            decode_mi(idx, m, k, a);
            mpfr_set(term, T[node * dcount + idx], RND);
            for (size_t d = 0; d < m; d++) mpfr_mul(term, term, factor[(d * 2 + c[d]) * (size_t)(k + 1) + a[d]], RND);
            mpfr_add(out, out, term, RND);
        }
    }
    for (size_t i = 0; i < m; i++) { mpfr_clear(h[i]); mpfr_clear(t[i]); }
    for (size_t i = 0; i < fcount; i++) mpfr_clear(factor[i]);
    mpfr_clears(phi, hp, term, tmp, (mpfr_ptr)0);
    free(lo); free(a); free(c); free(h); free(t); free(factor);
}

/* value-only tensor evaluation (default / spline), recursive over dims */
typedef struct { const MFun* f; mpfr_t* p; const int* ders; int64_t* s; size_t* w; int kernel; long prec; } MCtx;

static void eval_dim_m(const MCtx* ctx, size_t k, size_t base, mpfr_t* V, mpfr_t out) {
    const MFun* f = ctx->f;
    size_t wk = ctx->w[k];
    int64_t start = ctx->s[k];
    bool per = f->periodic && f->periodic[k];
    size_t ndist = per ? f->nk[k] - 1 : f->nk[k];

    mpfr_t* xs = malloc(sizeof(mpfr_t) * wk);
    mpfr_t* tmp = malloc(sizeof(mpfr_t) * wk);
    for (size_t t = 0; t < wk; t++) { mpfr_init2(xs[t], ctx->prec); mpfr_init2(tmp[t], ctx->prec); }
    for (size_t t = 0; t < wk; t++) {
        int64_t j = start + (int64_t)t;
        size_t r;
        if (per) {
            int64_t nd = (int64_t)ndist, rr = ((j % nd) + nd) % nd, q = (j - rr) / nd;
            r = (size_t)rr;
            mpfr_mul_si(xs[t], f->period[k], q, RND);   /* q * P */
            mpfr_add(xs[t], xs[t], f->grid[k][r], RND);
        } else { r = (size_t)j; mpfr_set(xs[t], f->grid[k][r], RND); }
        size_t off = base + r * f->stride[k];
        if (k + 1 == f->m) mpfr_set(tmp[t], V[off], RND);
        else eval_dim_m(ctx, k + 1, off, V, tmp[t]);
    }
    if (ctx->kernel == 1) {
        if (per) spline_eval_periodic_m(xs, wk, tmp, f->period[k], ctx->p[k], (size_t)ctx->ders[k], out, ctx->prec);
        else     spline_eval_m(xs, tmp, wk, ctx->p[k], (size_t)ctx->ders[k], out, ctx->prec);
    } else newton_eval_m(xs, tmp, wk, ctx->p[k], (size_t)ctx->ders[k], out, ctx->prec);
    for (size_t t = 0; t < wk; t++) { mpfr_clear(xs[t]); mpfr_clear(tmp[t]); }
    free(xs); free(tmp);
}

/* Evaluate one scalar component (selected by cpath) into the pre-init'd out. */
static bool eval_component_mpfr(MFun* f, const int* cpath, int vrank, size_t m,
                                mpfr_t* p, const int* ders, const int* orders,
                                int method, int Ksupplied, mpfr_t out, long prec) {
    if (method == METHOD_HERMITE || Ksupplied >= 1) {
        int K = Ksupplied, k = K > 1 ? K : 1;
        size_t dcount = ipow_sz((size_t)(k + 1), m);
        mpfr_t* basis = build_basis_m(k, prec);
        mpfr_t* T = build_T_m(f, m, K, k, dcount, cpath, vrank, prec);
        bool ok = (T != NULL);
        if (ok) {
            hermite_tensor_eval_m(f, T, dcount, m, k, basis, 2 * k + 2, p, ders, out, prec);
            for (size_t i = 0; i < f->total * dcount; i++) mpfr_clear(T[i]);
            free(T);
        }
        size_t nb = (size_t)2 * (k + 1) * (2 * k + 2);
        for (size_t i = 0; i < nb; i++) mpfr_clear(basis[i]);
        free(basis);
        return ok;
    }
    mpfr_t* V = malloc(sizeof(mpfr_t) * f->total);
    for (size_t i = 0; i < f->total; i++) mpfr_init2(V[i], prec);
    bool ok = true;
    for (size_t i = 0; i < f->total; i++)
        if (!value_component_m(f->entryAt[i], cpath, vrank, V[i])) { ok = false; break; }
    if (ok) {
        int64_t* s = malloc(sizeof(int64_t) * m);
        size_t* w = malloc(sizeof(size_t) * m);
        int kernel = (method == METHOD_SPLINE) ? 1 : 0;
        for (size_t k = 0; k < m; k++) {
            size_t nkk = f->nk[k];
            bool per = f->periodic && f->periodic[k];
            size_t ndist = per ? nkk - 1 : nkk;
            if (kernel == 1) { s[k] = 0; w[k] = ndist; continue; }
            size_t order;
            if (orders) order = (size_t)orders[k];
            else { order = nkk - 1; if (order > 3) order = 3; }
            if (order > ndist - 1) order = ndist - 1;
            size_t wk = order + 1;
            size_t i = bracket_m(f->grid[k], nkk, p[k]);
            size_t shift = (wk >= 2) ? (wk / 2 - 1) : 0;
            if (per) s[k] = (int64_t)i - (int64_t)shift;
            else { size_t sk = (i < shift) ? 0 : i - shift; if (sk > nkk - wk) sk = nkk - wk; s[k] = (int64_t)sk; }
            w[k] = wk;
        }
        MCtx ctx = { f, p, ders, s, w, kernel, prec };
        eval_dim_m(&ctx, 0, 0, V, out);
        free(s); free(w);
    }
    for (size_t i = 0; i < f->total; i++) mpfr_clear(V[i]);
    free(V);
    return ok;
}

/* Assemble a nested-List Expr from row-major mpfr components (moved, not freed). */
static Expr* assemble_array_m(mpfr_t* comps, const size_t* vshape, int vrank,
                              int level, size_t* idx) {
    if (level == vrank) return expr_new_mpfr_move(comps[(*idx)++]);
    size_t n = vshape[level];
    Expr** items = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) items[i] = assemble_array_m(comps, vshape, vrank, level + 1, idx);
    Expr* l = expr_new_function(expr_new_symbol("List"), items, n);
    free(items);
    return l;
}

/* --- public entry point ----------------------------------------------- */

Expr* interp_eval_mpfr(Expr* domain, Expr* table, size_t m,
                       Expr** call_args, const int* ders,
                       const int* orders, int method, int Ksupplied,
                       const bool* periodic, long bits) {
    long prec = bits > 53 ? bits : 53;

    mpfr_t* p = malloc(sizeof(mpfr_t) * m);
    for (size_t k = 0; k < m; k++) mpfr_init2(p[k], prec);
    bool okp = true;
    for (size_t k = 0; k < m; k++) if (!nm_set(call_args[k], p[k])) { okp = false; break; }
    if (!okp) { for (size_t k = 0; k < m; k++) mpfr_clear(p[k]); free(p); return NULL; }

    MFun f;
    if (!build_grid_m(domain, table, m, periodic, prec, &f)) {
        for (size_t k = 0; k < m; k++) mpfr_clear(p[k]);
        free(p);
        return NULL;
    }

    /* Reduce periodic coordinates into [x0, x0 + P). */
    mpfr_t u; mpfr_init2(u, prec);
    for (size_t k = 0; k < m; k++) {
        if (f.periodic[k]) {
            mpfr_sub(u, p[k], f.grid[k][0], RND);
            mpfr_fmod(u, u, f.period[k], RND);
            if (mpfr_sgn(u) < 0) mpfr_add(u, u, f.period[k], RND);
            mpfr_add(p[k], u, f.grid[k][0], RND);
        }
    }
    mpfr_clear(u);

    size_t vshape[16];
    int vrank;
    value_shape_m(m_entry_value(f.entryAt[0]), vshape, &vrank);
    size_t vtotal = 1;
    for (int dd = 0; dd < vrank; dd++) vtotal *= vshape[dd];

    mpfr_t* comps = malloc(sizeof(mpfr_t) * vtotal);
    for (size_t i = 0; i < vtotal; i++) mpfr_init2(comps[i], prec);
    int* cpath = malloc(sizeof(int) * (vrank ? (size_t)vrank : 1));
    bool ok = true;
    for (size_t ci = 0; ci < vtotal && ok; ci++) {
        size_t rem = ci;
        for (int dd = vrank; dd-- > 0; ) { cpath[dd] = (int)(rem % vshape[dd]); rem /= vshape[dd]; }
        ok = eval_component_mpfr(&f, cpath, vrank, m, p, ders, orders, method, Ksupplied, comps[ci], prec);
    }
    free(cpath);

    Expr* result = NULL;
    if (ok) {
        size_t idx = 0;
        result = assemble_array_m(comps, vshape, vrank, 0, &idx);   /* moves comps */
    } else {
        for (size_t i = 0; i < vtotal; i++) mpfr_clear(comps[i]);
    }
    free(comps);

    mfun_free(&f);
    for (size_t k = 0; k < m; k++) mpfr_clear(p[k]);
    free(p);
    return result;
}

#endif /* USE_MPFR */
