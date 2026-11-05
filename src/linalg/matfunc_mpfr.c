/*
 * matfunc_mpfr.c -- arbitrary-precision kernels for the matrix functions.
 *
 * Matrices are complex MPFR (MfMat: row-major re/im arrays).  A real input is
 * loaded with zero imaginary parts, and every algorithm below keeps an exactly
 * zero imaginary part zero, so a real answer comes back real.
 *
 *   mf_mpfr_expm   scaling and squaring with a truncated Taylor series: A is
 *                  scaled to ||2^-s A|| <= 2^-r with r ~ sqrt(prec) / 2, the
 *                  series is summed until a term drops below 2^-wprec of the
 *                  sum, and the result is squared s times.  The working
 *                  precision carries s + log2 n + 32 guard bits, so the answer
 *                  keeps the input's precision.
 *   mf_mpfr_sqrtm  the Denman--Beavers iteration in product form.
 *   mf_mpfr_logm   inverse scaling and squaring: Denman--Beavers square roots
 *                  until ||B - I|| <= 1/4, then log B = 2 atanh((B - I)(B + I)^-1)
 *                  by its odd series, times 2^k.
 *   mf_mpfr_powm   A^s = exp(s log A) (s = 1/2 by the square root).
 *
 * The square-root iteration converges to the PRINCIPAL root only when no
 * eigenvalue lies on the closed negative real axis.  The eigenvalues are
 * located by a machine-precision eigenvalue pass; when one is on (or near) the
 * cut, the matrix is rotated, A = e^(i theta) B with every eigenvalue of B
 * away from the cut, and log A = i theta I + log B (and sqrt A =
 * e^(i theta/2) sqrt B), which picks the same principal branch the reference
 * uses (Log[-1] = I Pi).  A zero eigenvalue returns MF_ZEROEIG.
 */

#include "matfunc_internal.h"
#include "lapack.h"
#include "numeric.h"
#include "sym_names.h"
#include "linalg.h"
#include "ndarray.h"
#include "la_dense.h"              /* la_has_complex */
#include "eval.h"
#include "numeric_complex.h"        /* numeric_mpfr_make_complex */
#include "schurdecomp_internal.h"  /* schur_mpfr_complex_core */

#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef USE_MPFR

/* ---- allocation ---------------------------------------------------------- */

static bool mat_init(MfMat* a, int n, mpfr_prec_t prec) {
    size_t nn = (size_t)n * (size_t)n;
    a->n = n; a->prec = prec;
    a->re = (mpfr_t*)malloc(sizeof(mpfr_t) * nn);
    a->im = (mpfr_t*)malloc(sizeof(mpfr_t) * nn);
    if (!a->re || !a->im) { free(a->re); free(a->im); a->re = a->im = NULL; return false; }
    for (size_t k = 0; k < nn; k++) {
        mpfr_init2(a->re[k], prec); mpfr_set_zero(a->re[k], 1);
        mpfr_init2(a->im[k], prec); mpfr_set_zero(a->im[k], 1);
    }
    return true;
}

void mf_mat_free(MfMat* a) {
    if (!a || !a->re) return;
    size_t nn = (size_t)a->n * (size_t)a->n;
    for (size_t k = 0; k < nn; k++) { mpfr_clear(a->re[k]); mpfr_clear(a->im[k]); }
    free(a->re); free(a->im);
    a->re = a->im = NULL;
}

static void mat_copy(MfMat* dst, const MfMat* src) {
    size_t nn = (size_t)src->n * (size_t)src->n;
    for (size_t k = 0; k < nn; k++) {
        mpfr_set(dst->re[k], src->re[k], MPFR_RNDN);
        mpfr_set(dst->im[k], src->im[k], MPFR_RNDN);
    }
}

static void mat_identity(MfMat* a, mpfr_t scale_re) {
    size_t nn = (size_t)a->n * (size_t)a->n;
    for (size_t k = 0; k < nn; k++) { mpfr_set_zero(a->re[k], 1); mpfr_set_zero(a->im[k], 1); }
    for (int i = 0; i < a->n; i++) {
        size_t d = (size_t)i * (size_t)a->n + (size_t)i;
        if (scale_re) mpfr_set(a->re[d], scale_re, MPFR_RNDN);
        else mpfr_set_ui(a->re[d], 1, MPFR_RNDN);
    }
}

#define IDX(a, i, j) ((size_t)(i) * (size_t)(a)->n + (size_t)(j))

bool mf_mat_load(const Expr* m, int n, mpfr_prec_t prec, MfMat* out) {
    if (!mat_init(out, n, prec)) return false;
    if (m->type != EXPR_FUNCTION || (int)m->data.function.arg_count != n) { mf_mat_free(out); return false; }
    for (int i = 0; i < n; i++) {
        const Expr* row = m->data.function.args[i];
        if (row->type != EXPR_FUNCTION || (int)row->data.function.arg_count != n) { mf_mat_free(out); return false; }
        for (int j = 0; j < n; j++) {
            bool inexact = false;
            if (!get_approx_mpfr(row->data.function.args[j], out->re[IDX(out, i, j)],
                                 out->im[IDX(out, i, j)], &inexact)) {
                mf_mat_free(out);
                return false;
            }
        }
    }
    return true;
}

Expr* mf_mat_to_expr(const MfMat* a, mpfr_prec_t out_bits, bool force_real) {
    int n = a->n;
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    Expr** cells = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    if (!rows || !cells) { free(rows); free(cells); return NULL; }
    bool any_im = false;
    if (!force_real) {
        size_t nn = (size_t)n * (size_t)n;
        for (size_t k = 0; k < nn && !any_im; k++) if (!mpfr_zero_p(a->im[k])) any_im = true;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mpfr_t r;
            mpfr_init2(r, out_bits);
            mpfr_set(r, a->re[IDX(a, i, j)], MPFR_RNDN);
            Expr* re = expr_new_mpfr_move(r);
            if (any_im) {
                mpfr_t im;
                mpfr_init2(im, out_bits);
                mpfr_set(im, a->im[IDX(a, i, j)], MPFR_RNDN);
                Expr* args[2] = { re, expr_new_mpfr_move(im) };
                cells[j] = expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
            } else {
                cells[j] = re;
            }
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)n);
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows); free(cells);
    return r;
}

/* ---- arithmetic ---------------------------------------------------------- */

/* C = A B (C distinct from A, B). */
static void mat_mul(MfMat* C, const MfMat* A, const MfMat* B) {
    int n = A->n;
    mpfr_t sr, si, t;
    mpfr_inits2(C->prec, sr, si, t, (mpfr_ptr)0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            mpfr_set_zero(sr, 1); mpfr_set_zero(si, 1);
            for (int k = 0; k < n; k++) {
                size_t a = IDX(A, i, k), b = IDX(B, k, j);
                mpfr_mul(t, A->re[a], B->re[b], MPFR_RNDN); mpfr_add(sr, sr, t, MPFR_RNDN);
                if (!mpfr_zero_p(A->im[a]) && !mpfr_zero_p(B->im[b])) {
                    mpfr_mul(t, A->im[a], B->im[b], MPFR_RNDN); mpfr_sub(sr, sr, t, MPFR_RNDN);
                }
                if (!mpfr_zero_p(B->im[b])) {
                    mpfr_mul(t, A->re[a], B->im[b], MPFR_RNDN); mpfr_add(si, si, t, MPFR_RNDN);
                }
                if (!mpfr_zero_p(A->im[a])) {
                    mpfr_mul(t, A->im[a], B->re[b], MPFR_RNDN); mpfr_add(si, si, t, MPFR_RNDN);
                }
            }
            mpfr_set(C->re[IDX(C, i, j)], sr, MPFR_RNDN);
            mpfr_set(C->im[IDX(C, i, j)], si, MPFR_RNDN);
        }
    mpfr_clears(sr, si, t, (mpfr_ptr)0);
}

/* || A ||_1 as a double (an estimate is all the callers need). */
static double mat_norm1(const MfMat* A) {
    double best = 0.0;
    for (int j = 0; j < A->n; j++) {
        double s = 0.0;
        for (int i = 0; i < A->n; i++)
            s += hypot(mpfr_get_d(A->re[IDX(A, i, j)], MPFR_RNDN),
                       mpfr_get_d(A->im[IDX(A, i, j)], MPFR_RNDN));
        if (s > best) best = s;
    }
    return best;
}

/* log2 || A ||_1 without overflowing a double (for huge MPFR exponents). */
static double mat_log2norm(const MfMat* A) {
    double best = -INFINITY;
    for (int j = 0; j < A->n; j++) {
        double s = 0.0;
        long emax = LONG_MIN;
        for (int i = 0; i < A->n; i++) {
            if (!mpfr_zero_p(A->re[IDX(A, i, j)])) { long e = mpfr_get_exp(A->re[IDX(A, i, j)]); if (e > emax) emax = e; }
            if (!mpfr_zero_p(A->im[IDX(A, i, j)])) { long e = mpfr_get_exp(A->im[IDX(A, i, j)]); if (e > emax) emax = e; }
        }
        if (emax == LONG_MIN) continue;
        for (int i = 0; i < A->n; i++) {
            long e1 = 0, e2 = 0;
            double a = mpfr_get_d_2exp(&e1, A->re[IDX(A, i, j)], MPFR_RNDN);
            double b = mpfr_get_d_2exp(&e2, A->im[IDX(A, i, j)], MPFR_RNDN);
            s += fabs(ldexp(a, (int)(e1 - emax))) + fabs(ldexp(b, (int)(e2 - emax)));
        }
        double l = log2(s) + (double)emax;
        if (l > best) best = l;
    }
    return best;
}

/* Inverse by Gauss--Jordan with partial pivoting; false when singular. */
static bool mat_inverse(MfMat* X, const MfMat* A) {
    int n = A->n;
    MfMat W;
    if (!mat_init(&W, n, A->prec)) return false;
    mat_copy(&W, A);
    mat_identity(X, NULL);
    mpfr_t tr, ti, ar, ai, br, bi, d, mag, best;
    mpfr_inits2(A->prec, tr, ti, ar, ai, br, bi, d, mag, best, (mpfr_ptr)0);
    bool ok = true;
    for (int k = 0; k < n && ok; k++) {
        int piv = -1;
        mpfr_set_si(best, -1, MPFR_RNDN);
        for (int i = k; i < n; i++) {
            mpfr_hypot(mag, W.re[IDX(&W, i, k)], W.im[IDX(&W, i, k)], MPFR_RNDN);
            if (mpfr_cmp(mag, best) > 0) { mpfr_set(best, mag, MPFR_RNDN); piv = i; }
        }
        if (piv < 0 || mpfr_zero_p(best)) { ok = false; break; }
        if (piv != k)
            for (int j = 0; j < n; j++) {
                mpfr_swap(W.re[IDX(&W, k, j)], W.re[IDX(&W, piv, j)]);
                mpfr_swap(W.im[IDX(&W, k, j)], W.im[IDX(&W, piv, j)]);
                mpfr_swap(X->re[IDX(X, k, j)], X->re[IDX(X, piv, j)]);
                mpfr_swap(X->im[IDX(X, k, j)], X->im[IDX(X, piv, j)]);
            }
        /* (ar, ai) = 1 / W[k][k] */
        mpfr_sqr(d, W.re[IDX(&W, k, k)], MPFR_RNDN);
        mpfr_sqr(tr, W.im[IDX(&W, k, k)], MPFR_RNDN);
        mpfr_add(d, d, tr, MPFR_RNDN);
        mpfr_div(ar, W.re[IDX(&W, k, k)], d, MPFR_RNDN);
        mpfr_div(ai, W.im[IDX(&W, k, k)], d, MPFR_RNDN);
        mpfr_neg(ai, ai, MPFR_RNDN);
        /* scale row k */
        for (int pass = 0; pass < 2; pass++) {
            MfMat* M = pass == 0 ? &W : X;
            for (int j = 0; j < n; j++) {
                size_t q = IDX(M, k, j);
                mpfr_mul(tr, M->re[q], ar, MPFR_RNDN);
                mpfr_mul(ti, M->im[q], ai, MPFR_RNDN);
                mpfr_sub(br, tr, ti, MPFR_RNDN);
                mpfr_mul(tr, M->re[q], ai, MPFR_RNDN);
                mpfr_mul(ti, M->im[q], ar, MPFR_RNDN);
                mpfr_add(bi, tr, ti, MPFR_RNDN);
                mpfr_set(M->re[q], br, MPFR_RNDN);
                mpfr_set(M->im[q], bi, MPFR_RNDN);
            }
        }
        /* eliminate column k from every other row */
        for (int i = 0; i < n; i++) {
            if (i == k) continue;
            mpfr_set(br, W.re[IDX(&W, i, k)], MPFR_RNDN);
            mpfr_set(bi, W.im[IDX(&W, i, k)], MPFR_RNDN);
            if (mpfr_zero_p(br) && mpfr_zero_p(bi)) continue;
            for (int pass = 0; pass < 2; pass++) {
                MfMat* M = pass == 0 ? &W : X;
                for (int j = 0; j < n; j++) {
                    size_t qk = IDX(M, k, j), qi = IDX(M, i, j);
                    mpfr_mul(tr, br, M->re[qk], MPFR_RNDN);
                    mpfr_mul(ti, bi, M->im[qk], MPFR_RNDN);
                    mpfr_sub(tr, tr, ti, MPFR_RNDN);
                    mpfr_sub(M->re[qi], M->re[qi], tr, MPFR_RNDN);
                    mpfr_mul(tr, br, M->im[qk], MPFR_RNDN);
                    mpfr_mul(ti, bi, M->re[qk], MPFR_RNDN);
                    mpfr_add(tr, tr, ti, MPFR_RNDN);
                    mpfr_sub(M->im[qi], M->im[qi], tr, MPFR_RNDN);
                }
            }
        }
    }
    mpfr_clears(tr, ti, ar, ai, br, bi, d, mag, best, (mpfr_ptr)0);
    mf_mat_free(&W);
    return ok;
}

/* Raise the working precision of a matrix (new storage). */
static bool mat_reprec(MfMat* dst, const MfMat* src, mpfr_prec_t prec) {
    if (!mat_init(dst, src->n, prec)) return false;
    mat_copy(dst, src);
    return true;
}

/* Round `w` back into `a` at a's own precision. */
static void mat_store(MfMat* a, const MfMat* w) { mat_copy(a, w); }

static mpfr_prec_t guard_prec(const MfMat* a, long extra) {
    long lg = 0;
    while ((1L << lg) < a->n + 1) lg++;
    return a->prec + 32 + 2 * lg + extra;
}

/* ---- exp ----------------------------------------------------------------- */

static int expm_work(MfMat* W) {
    int n = W->n;
    mpfr_prec_t wp = W->prec;
    double l2 = mat_log2norm(W);
    if (l2 == -INFINITY) {             /* exp(0) = I */
        mat_identity(W, NULL);
        return MF_OK;
    }
    if (!isfinite(l2) || l2 > 1e6) return MF_FAIL;
    long r = (long)ceil(sqrt((double)wp) / 2.0);
    long s = (long)ceil(l2) + r;
    if (s < 0) s = 0;
    /* B = 2^-s A */
    size_t nn = (size_t)n * (size_t)n;
    for (size_t k = 0; k < nn; k++) {
        mpfr_div_2si(W->re[k], W->re[k], s, MPFR_RNDN);
        mpfr_div_2si(W->im[k], W->im[k], s, MPFR_RNDN);
    }
    MfMat S, term, tmp;
    if (!mat_init(&S, n, wp)) return MF_FAIL;
    if (!mat_init(&term, n, wp)) { mf_mat_free(&S); return MF_FAIL; }
    if (!mat_init(&tmp, n, wp)) { mf_mat_free(&S); mf_mat_free(&term); return MF_FAIL; }
    mat_identity(&S, NULL);
    mat_identity(&term, NULL);
    for (long k = 1; k < 100000; k++) {
        mat_mul(&tmp, &term, W);
        for (size_t q = 0; q < nn; q++) {
            mpfr_div_si(term.re[q], tmp.re[q], k, MPFR_RNDN);
            mpfr_div_si(term.im[q], tmp.im[q], k, MPFR_RNDN);
            mpfr_add(S.re[q], S.re[q], term.re[q], MPFR_RNDN);
            mpfr_add(S.im[q], S.im[q], term.im[q], MPFR_RNDN);
        }
        double tl = mat_log2norm(&term), sl = mat_log2norm(&S);
        if (tl == -INFINITY || tl < sl - (double)wp - 4.0) break;
    }
    for (long q = 0; q < s; q++) {
        mat_mul(&tmp, &S, &S);
        mat_copy(&S, &tmp);
    }
    mat_copy(W, &S);
    mf_mat_free(&S); mf_mat_free(&term); mf_mat_free(&tmp);
    return MF_OK;
}

int mf_mpfr_expm(MfMat* a) {
    double l2 = mat_log2norm(a);
    long s_est = (l2 > 0.0 && isfinite(l2)) ? (long)ceil(l2) : 0;
    long r = (long)ceil(sqrt((double)a->prec) / 2.0);
    MfMat W;
    if (!mat_reprec(&W, a, guard_prec(a, s_est + r))) return MF_FAIL;
    int st = expm_work(&W);
    if (st == MF_OK) mat_store(a, &W);
    mf_mat_free(&W);
    return st;
}

/* ---- eigenvalue location (for the branch cut) ----------------------------- */

/* theta (radians) so that every eigenvalue of e^(-i theta) A lies off the
 * closed negative real axis; 0 when none is near the cut.  MF_ZEROEIG for a
 * (numerically) zero eigenvalue, MF_FAIL when no rotation separates them or
 * no eigenvalue solver is available. */
static int branch_rotation(const MfMat* a, double* theta) {
    *theta = 0.0;
    int n = a->n;
    size_t nn = (size_t)n * (size_t)n;
    double* Z = (double*)malloc(sizeof(double) * 2 * nn);
    double* w = (double*)malloc(sizeof(double) * 2 * (size_t)n);
    if (!Z || !w) { free(Z); free(w); return MF_FAIL; }
    double nrm = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            size_t c = (size_t)i + (size_t)j * (size_t)n;
            Z[2 * c]     = mpfr_get_d(a->re[IDX(a, i, j)], MPFR_RNDN);
            Z[2 * c + 1] = mpfr_get_d(a->im[IDX(a, i, j)], MPFR_RNDN);
            nrm = fmax(nrm, hypot(Z[2 * c], Z[2 * c + 1]));
        }
    int st = MF_OK;
    if (!mathilda_lapack_probe() || mat_lapack_zgeev_values(n, Z, n, w) != 0) {
        /* No eigenvalue pass: assume no rotation is needed (Denman--Beavers
         * then simply fails to converge for an eigenvalue on the cut). */
        free(Z); free(w);
        return MF_OK;
    }
    double min_arg = M_PI;
    bool near_cut = false;
    for (int i = 0; i < n; i++) {
        double lr = w[2 * i], li = w[2 * i + 1];
        double la = hypot(lr, li);
        if (la <= 1e-13 * (nrm > 0.0 ? nrm : 1.0) * n) { st = MF_ZEROEIG; break; }
        double arg = (fabs(li) <= 1e-12 * la && lr < 0.0) ? M_PI : atan2(li, lr);
        if (arg < min_arg) min_arg = arg;
        if (fabs(arg) > M_PI - 0.25) near_cut = true;
    }
    if (st == MF_OK && near_cut) {
        double margin = 0.5 * (M_PI + min_arg);
        if (margin < 1e-6) st = MF_FAIL;
        else *theta = margin;
    }
    free(Z); free(w);
    return st;
}

/* W *= e^(-i theta) (in place). */
static void rotate(MfMat* W, double theta, int sign) {
    if (theta == 0.0) return;
    mpfr_t c, s, t1, t2, nr, pi;
    mpfr_inits2(W->prec, c, s, t1, t2, nr, pi, (mpfr_ptr)0);
    /* theta is a double; use it exactly (any angle away from the cut works,
     * so the rotation itself need not be a "nice" number). */
    mpfr_set_d(t1, theta * (double)sign, MPFR_RNDN);
    mpfr_sin_cos(s, c, t1, MPFR_RNDN);
    mpfr_neg(s, s, MPFR_RNDN);        /* e^(-i sign theta) = c - i s' */
    size_t nn = (size_t)W->n * (size_t)W->n;
    for (size_t k = 0; k < nn; k++) {
        mpfr_mul(t1, W->re[k], c, MPFR_RNDN);
        mpfr_mul(t2, W->im[k], s, MPFR_RNDN);
        mpfr_sub(nr, t1, t2, MPFR_RNDN);
        mpfr_mul(t1, W->re[k], s, MPFR_RNDN);
        mpfr_mul(t2, W->im[k], c, MPFR_RNDN);
        mpfr_add(W->im[k], t1, t2, MPFR_RNDN);
        mpfr_set(W->re[k], nr, MPFR_RNDN);
    }
    mpfr_clears(c, s, t1, t2, nr, pi, (mpfr_ptr)0);
}

/* ---- square root (Denman--Beavers, product form) -------------------------- */

static int sqrtm_db(MfMat* W) {
    int n = W->n;
    size_t nn = (size_t)n * (size_t)n;
    MfMat M, Y, Mi, T;
    mpfr_prec_t wp = W->prec;
    if (!mat_init(&M, n, wp)) return MF_FAIL;
    if (!mat_init(&Y, n, wp)) { mf_mat_free(&M); return MF_FAIL; }
    if (!mat_init(&Mi, n, wp)) { mf_mat_free(&M); mf_mat_free(&Y); return MF_FAIL; }
    if (!mat_init(&T, n, wp)) { mf_mat_free(&M); mf_mat_free(&Y); mf_mat_free(&Mi); return MF_FAIL; }
    mat_copy(&M, W);
    mat_copy(&Y, W);
    int st = MF_FAIL;
    for (int it = 0; it < 200; it++) {
        if (!mat_inverse(&Mi, &M)) { st = MF_ZEROEIG; break; }
        /* Y <- Y (I + Mi) / 2 */
        for (int i = 0; i < n; i++) mpfr_add_ui(Mi.re[IDX(&Mi, i, i)], Mi.re[IDX(&Mi, i, i)], 1, MPFR_RNDN);
        mat_mul(&T, &Y, &Mi);
        for (size_t k = 0; k < nn; k++) {
            mpfr_div_2ui(Y.re[k], T.re[k], 1, MPFR_RNDN);
            mpfr_div_2ui(Y.im[k], T.im[k], 1, MPFR_RNDN);
        }
        /* M <- (I + (M + Mi) / 2) / 2  ==  M/4 + Mi'/4 + I/4, Mi' = Mi - I + ... */
        for (int i = 0; i < n; i++) mpfr_sub_ui(Mi.re[IDX(&Mi, i, i)], Mi.re[IDX(&Mi, i, i)], 1, MPFR_RNDN);
        for (size_t k = 0; k < nn; k++) {
            mpfr_add(M.re[k], M.re[k], Mi.re[k], MPFR_RNDN);
            mpfr_add(M.im[k], M.im[k], Mi.im[k], MPFR_RNDN);
            mpfr_div_2ui(M.re[k], M.re[k], 2, MPFR_RNDN);
            mpfr_div_2ui(M.im[k], M.im[k], 2, MPFR_RNDN);
        }
        mpfr_t half;
        mpfr_init2(half, wp);
        mpfr_set_d(half, 0.5, MPFR_RNDN);
        for (int i = 0; i < n; i++) mpfr_add(M.re[IDX(&M, i, i)], M.re[IDX(&M, i, i)], half, MPFR_RNDN);
        mpfr_clear(half);
        /* converged when ||M - I|| is below the working precision */
        mat_copy(&T, &M);
        for (int i = 0; i < n; i++) mpfr_sub_ui(T.re[IDX(&T, i, i)], T.re[IDX(&T, i, i)], 1, MPFR_RNDN);
        double l2 = mat_log2norm(&T);
        if (!isfinite(l2) && l2 > 0) break;
        if (l2 == -INFINITY || l2 < -(double)wp + 8.0) { st = MF_OK; break; }
    }
    if (st == MF_OK) mat_copy(W, &Y);
    mf_mat_free(&M); mf_mat_free(&Y); mf_mat_free(&Mi); mf_mat_free(&T);
    return st;
}

int mf_mpfr_sqrtm(MfMat* a) {
    double theta = 0.0;
    int st = branch_rotation(a, &theta);
    if (st != MF_OK) return st;
    MfMat W;
    if (!mat_reprec(&W, a, guard_prec(a, 16))) return MF_FAIL;
    rotate(&W, theta, 1);
    st = sqrtm_db(&W);
    if (st == MF_OK) {
        rotate(&W, theta / 2.0, -1);   /* sqrt A = e^(i theta/2) sqrt B */
        mat_store(a, &W);
    }
    mf_mat_free(&W);
    return st;
}

/* ---- log ------------------------------------------------------------------ */

/* log of W (no eigenvalue on the cut) in place, by inverse scaling and squaring. */
static int logm_work(MfMat* W) {
    int n = W->n;
    size_t nn = (size_t)n * (size_t)n;
    mpfr_prec_t wp = W->prec;
    MfMat X;
    if (!mat_init(&X, n, wp)) return MF_FAIL;
    int k = 0, st = MF_OK;
    for (;;) {
        mat_copy(&X, W);
        for (int i = 0; i < n; i++) mpfr_sub_ui(X.re[IDX(&X, i, i)], X.re[IDX(&X, i, i)], 1, MPFR_RNDN);
        if (mat_norm1(&X) <= 0.25) break;
        if (k > 200) { st = MF_FAIL; break; }
        st = sqrtm_db(W);
        if (st != MF_OK) break;
        k++;
    }
    if (st == MF_OK) {
        /* Z = (W - I)(W + I)^-1 ; log W = 2 sum_{j odd} Z^j / j */
        MfMat P, Z, Zsq, S, T;
        mat_init(&P, n, wp); mat_init(&Z, n, wp); mat_init(&Zsq, n, wp);
        mat_init(&S, n, wp); mat_init(&T, n, wp);
        mat_copy(&P, W);
        for (int i = 0; i < n; i++) mpfr_add_ui(P.re[IDX(&P, i, i)], P.re[IDX(&P, i, i)], 1, MPFR_RNDN);
        if (!mat_inverse(&T, &P)) st = MF_FAIL;
        if (st == MF_OK) {
            mat_mul(&Z, &X, &T);           /* X = W - I */
            mat_mul(&Zsq, &Z, &Z);
            mat_copy(&S, &Z);
            mat_copy(&P, &Z);              /* P = Z^j */
            for (long j = 3; j < 1000000; j += 2) {
                mat_mul(&T, &P, &Zsq);
                mat_copy(&P, &T);
                double tl = mat_log2norm(&P) - log2((double)j);
                for (size_t q = 0; q < nn; q++) {
                    mpfr_div_si(T.re[q], P.re[q], j, MPFR_RNDN);
                    mpfr_div_si(T.im[q], P.im[q], j, MPFR_RNDN);
                    mpfr_add(S.re[q], S.re[q], T.re[q], MPFR_RNDN);
                    mpfr_add(S.im[q], S.im[q], T.im[q], MPFR_RNDN);
                }
                double sl = mat_log2norm(&S);
                if (tl == -INFINITY || tl < sl - (double)wp - 4.0) break;
            }
            for (size_t q = 0; q < nn; q++) {
                mpfr_mul_2si(W->re[q], S.re[q], k + 1, MPFR_RNDN);
                mpfr_mul_2si(W->im[q], S.im[q], k + 1, MPFR_RNDN);
            }
        }
        mf_mat_free(&P); mf_mat_free(&Z); mf_mat_free(&Zsq); mf_mat_free(&S); mf_mat_free(&T);
    }
    mf_mat_free(&X);
    return st;
}

int mf_mpfr_logm(MfMat* a) {
    double theta = 0.0;
    int st = branch_rotation(a, &theta);
    if (st != MF_OK) return st;
    MfMat W;
    if (!mat_reprec(&W, a, guard_prec(a, 32))) return MF_FAIL;
    rotate(&W, theta, 1);
    st = logm_work(&W);
    if (st == MF_OK) {
        if (theta != 0.0) {              /* log A = i theta I + log B */
            mpfr_t t;
            mpfr_init2(t, W.prec);
            mpfr_set_d(t, theta, MPFR_RNDN);
            for (int i = 0; i < W.n; i++)
                mpfr_add(W.im[IDX(&W, i, i)], W.im[IDX(&W, i, i)], t, MPFR_RNDN);
            mpfr_clear(t);
        }
        mat_store(a, &W);
    }
    mf_mat_free(&W);
    return st;
}

int mf_mpfr_powm(MfMat* a, const mpfr_t s_re, const mpfr_t s_im) {
    if (mpfr_zero_p(s_im) && mpfr_cmp_d(s_re, 0.5) == 0) return mf_mpfr_sqrtm(a);
    double theta = 0.0;
    int st = branch_rotation(a, &theta);
    if (st != MF_OK) return st;
    double l2 = mat_log2norm(a);
    long extra = 64 + (l2 > 0.0 && isfinite(l2) ? (long)ceil(log2(l2 + 1.0)) : 0);
    MfMat W;
    if (!mat_reprec(&W, a, guard_prec(a, extra))) return MF_FAIL;
    rotate(&W, theta, 1);
    st = logm_work(&W);
    if (st == MF_OK) {
        if (theta != 0.0) {
            mpfr_t t;
            mpfr_init2(t, W.prec);
            mpfr_set_d(t, theta, MPFR_RNDN);
            for (int i = 0; i < W.n; i++)
                mpfr_add(W.im[IDX(&W, i, i)], W.im[IDX(&W, i, i)], t, MPFR_RNDN);
            mpfr_clear(t);
        }
        /* W <- s W */
        size_t nn = (size_t)W.n * (size_t)W.n;
        mpfr_t t1, t2, nr;
        mpfr_inits2(W.prec, t1, t2, nr, (mpfr_ptr)0);
        for (size_t q = 0; q < nn; q++) {
            mpfr_mul(t1, W.re[q], s_re, MPFR_RNDN);
            mpfr_mul(t2, W.im[q], s_im, MPFR_RNDN);
            mpfr_sub(nr, t1, t2, MPFR_RNDN);
            mpfr_mul(t1, W.re[q], s_im, MPFR_RNDN);
            mpfr_mul(t2, W.im[q], s_re, MPFR_RNDN);
            mpfr_add(W.im[q], t1, t2, MPFR_RNDN);
            mpfr_set(W.re[q], nr, MPFR_RNDN);
        }
        mpfr_clears(t1, t2, nr, (mpfr_ptr)0);
        st = expm_work(&W);
        if (st == MF_OK) mat_store(a, &W);
    }
    mf_mat_free(&W);
    return st;
}

/* f(A) of an arbitrary-precision matrix by the Parlett recurrence on its
 * complex Schur form A = Q T Q^H (schur_mpfr_complex_core):
 *     F_ii = f(T_ii),
 *     F_ij = [T_ij (F_jj - F_ii) + sum_{i<k<j} (T_ik F_kj - F_ik T_kj)]
 *            / (T_jj - T_ii),
 * then F(A) = Q F Q^H.  The work runs at twice the input precision plus a
 * guard, so the divided differences of eigenvalues as close as 2^-bits keep
 * the input's digits; eigenvalues closer than that (a repeated eigenvalue, a
 * derogatory matrix) decline to the exact path.  It replaced that exact path
 * for every MPFR matrix: rationalising a 6 x 6 30-digit matrix and taking
 * MatrixFunction[Sin, .] symbolically ran past any time limit. */
static bool mpfr_scalar_f(const Expr* f, mpfr_t zr, mpfr_t zi, mpfr_t fr, mpfr_t fi)
{
    Expr* z = mpfr_zero_p(zi) ? expr_new_mpfr_copy(zr) : numeric_mpfr_make_complex(zr, zi);
    Expr* a[1] = { z };
    Expr* val = eval_and_free(expr_new_function(expr_copy((Expr*)f), a, 1));
    bool inexact = false;
    mpfr_set_zero(fi, 1);
    bool ok = get_approx_mpfr(val, fr, fi, &inexact)
              && mpfr_number_p(fr) && mpfr_number_p(fi);
    expr_free(val);
    return ok;
}

Expr* mf_mpfr_function(const Expr* m, int n, long bits_in, const Expr* f)
{
    mpfr_prec_t bits = (mpfr_prec_t)bits_in;
    mpfr_prec_t wp = 2 * bits + 32;
    mpfr_t *Tr = NULL, *Ti = NULL, *Qr = NULL, *Qi = NULL;
    if (!schur_mpfr_complex_core(m, n, wp, &Tr, &Ti, &Qr, &Qi)) return NULL;
    size_t nn = (size_t)n * (size_t)n;
#define AT(i, j) ((size_t)(i) * (size_t)n + (size_t)(j))
    mpfr_t* Fr = (mpfr_t*)malloc(sizeof(mpfr_t) * nn);
    mpfr_t* Fi = (mpfr_t*)malloc(sizeof(mpfr_t) * nn);
    for (size_t k = 0; k < nn; k++) {
        mpfr_init2(Fr[k], wp); mpfr_set_zero(Fr[k], 1);
        mpfr_init2(Fi[k], wp); mpfr_set_zero(Fi[k], 1);
    }
    mpfr_t sr, si, tr, ti, dr, di, den, sep, scale;
    mpfr_inits2(wp, sr, si, tr, ti, dr, di, den, sep, scale, (mpfr_ptr)0);
    bool ok = true;
    /* the eigenvalue scale, for the separation test */
    mpfr_set_ui(scale, 1, MPFR_RNDN);
    for (int i = 0; i < n; i++) {
        mpfr_hypot(tr, Tr[AT(i, i)], Ti[AT(i, i)], MPFR_RNDN);
        if (mpfr_cmp(tr, scale) > 0) mpfr_set(scale, tr, MPFR_RNDN);
    }
    mpfr_mul_2si(sep, scale, -(long)bits, MPFR_RNDN);
    for (int i = 0; i < n && ok; i++)
        ok = mpfr_scalar_f(f, Tr[AT(i, i)], Ti[AT(i, i)], Fr[AT(i, i)], Fi[AT(i, i)]);
    for (int d = 1; d < n && ok; d++) {
        for (int i = 0; i + d < n && ok; i++) {
            int j = i + d;
            /* s = T_ij (F_jj - F_ii) */
            mpfr_sub(dr, Fr[AT(j, j)], Fr[AT(i, i)], MPFR_RNDN);
            mpfr_sub(di, Fi[AT(j, j)], Fi[AT(i, i)], MPFR_RNDN);
            mpfr_mul(sr, Tr[AT(i, j)], dr, MPFR_RNDN);
            mpfr_mul(tr, Ti[AT(i, j)], di, MPFR_RNDN); mpfr_sub(sr, sr, tr, MPFR_RNDN);
            mpfr_mul(si, Tr[AT(i, j)], di, MPFR_RNDN);
            mpfr_mul(tr, Ti[AT(i, j)], dr, MPFR_RNDN); mpfr_add(si, si, tr, MPFR_RNDN);
            for (int k = i + 1; k < j; k++) {
                /* + T_ik F_kj - F_ik T_kj   (F T = T F, entry (i, j)) */
                mpfr_mul(tr, Fr[AT(i, k)], Tr[AT(k, j)], MPFR_RNDN); mpfr_sub(sr, sr, tr, MPFR_RNDN);
                mpfr_mul(tr, Fi[AT(i, k)], Ti[AT(k, j)], MPFR_RNDN); mpfr_add(sr, sr, tr, MPFR_RNDN);
                mpfr_mul(tr, Fr[AT(i, k)], Ti[AT(k, j)], MPFR_RNDN); mpfr_sub(si, si, tr, MPFR_RNDN);
                mpfr_mul(tr, Fi[AT(i, k)], Tr[AT(k, j)], MPFR_RNDN); mpfr_sub(si, si, tr, MPFR_RNDN);
                mpfr_mul(tr, Tr[AT(i, k)], Fr[AT(k, j)], MPFR_RNDN); mpfr_add(sr, sr, tr, MPFR_RNDN);
                mpfr_mul(tr, Ti[AT(i, k)], Fi[AT(k, j)], MPFR_RNDN); mpfr_sub(sr, sr, tr, MPFR_RNDN);
                mpfr_mul(tr, Tr[AT(i, k)], Fi[AT(k, j)], MPFR_RNDN); mpfr_add(si, si, tr, MPFR_RNDN);
                mpfr_mul(tr, Ti[AT(i, k)], Fr[AT(k, j)], MPFR_RNDN); mpfr_add(si, si, tr, MPFR_RNDN);
            }
            /* / (T_jj - T_ii) */
            mpfr_sub(dr, Tr[AT(j, j)], Tr[AT(i, i)], MPFR_RNDN);
            mpfr_sub(di, Ti[AT(j, j)], Ti[AT(i, i)], MPFR_RNDN);
            mpfr_hypot(den, dr, di, MPFR_RNDN);
            if (mpfr_cmp(den, sep) <= 0) {
                /* a repeated (or numerically repeated) eigenvalue: the
                 * divided difference is a derivative the recurrence cannot
                 * form -- unless the block is structurally decoupled (T_ij
                 * and every T_ik, T_kj between them exactly 0, as for a
                 * diagonal matrix), when F_ij is 0 */
                bool decoupled = mpfr_zero_p(Tr[AT(i, j)]) && mpfr_zero_p(Ti[AT(i, j)]);
                for (int k = i + 1; k < j && decoupled; k++)
                    decoupled = mpfr_zero_p(Tr[AT(i, k)]) && mpfr_zero_p(Ti[AT(i, k)])
                             && mpfr_zero_p(Tr[AT(k, j)]) && mpfr_zero_p(Ti[AT(k, j)]);
                if (decoupled) {
                    mpfr_set_zero(Fr[AT(i, j)], 1); mpfr_set_zero(Fi[AT(i, j)], 1);
                    continue;
                }
                ok = false; break;
            }
            mpfr_sqr(den, den, MPFR_RNDN);
            /* (s)(conj d) / |d|^2 */
            mpfr_mul(tr, sr, dr, MPFR_RNDN); mpfr_mul(ti, si, di, MPFR_RNDN);
            mpfr_add(tr, tr, ti, MPFR_RNDN);
            mpfr_div(Fr[AT(i, j)], tr, den, MPFR_RNDN);
            mpfr_mul(tr, si, dr, MPFR_RNDN); mpfr_mul(ti, sr, di, MPFR_RNDN);
            mpfr_sub(tr, tr, ti, MPFR_RNDN);
            mpfr_div(Fi[AT(i, j)], tr, den, MPFR_RNDN);
        }
    }
    Expr* out = NULL;
    if (ok) {
        /* R = Q F Q^H, into the T arrays (no longer needed) */
        MfMat W;
        if (mat_init(&W, n, bits)) {
            mpfr_t ar, ai;
            mpfr_inits2(wp, ar, ai, (mpfr_ptr)0);
            /* G = Q F  (into T) */
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    mpfr_set_zero(ar, 1); mpfr_set_zero(ai, 1);
                    for (int k = 0; k <= j; k++) {
                        mpfr_mul(tr, Qr[AT(i, k)], Fr[AT(k, j)], MPFR_RNDN); mpfr_add(ar, ar, tr, MPFR_RNDN);
                        mpfr_mul(tr, Qi[AT(i, k)], Fi[AT(k, j)], MPFR_RNDN); mpfr_sub(ar, ar, tr, MPFR_RNDN);
                        mpfr_mul(tr, Qr[AT(i, k)], Fi[AT(k, j)], MPFR_RNDN); mpfr_add(ai, ai, tr, MPFR_RNDN);
                        mpfr_mul(tr, Qi[AT(i, k)], Fr[AT(k, j)], MPFR_RNDN); mpfr_add(ai, ai, tr, MPFR_RNDN);
                    }
                    mpfr_set(Tr[AT(i, j)], ar, MPFR_RNDN);
                    mpfr_set(Ti[AT(i, j)], ai, MPFR_RNDN);
                }
            /* R = G Q^H */
            bool real_in = !la_has_complex(m);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) {
                    mpfr_set_zero(ar, 1); mpfr_set_zero(ai, 1);
                    for (int k = 0; k < n; k++) {
                        /* G_ik conj(Q_jk) */
                        mpfr_mul(tr, Tr[AT(i, k)], Qr[AT(j, k)], MPFR_RNDN); mpfr_add(ar, ar, tr, MPFR_RNDN);
                        mpfr_mul(tr, Ti[AT(i, k)], Qi[AT(j, k)], MPFR_RNDN); mpfr_add(ar, ar, tr, MPFR_RNDN);
                        mpfr_mul(tr, Ti[AT(i, k)], Qr[AT(j, k)], MPFR_RNDN); mpfr_add(ai, ai, tr, MPFR_RNDN);
                        mpfr_mul(tr, Tr[AT(i, k)], Qi[AT(j, k)], MPFR_RNDN); mpfr_sub(ai, ai, tr, MPFR_RNDN);
                    }
                    mpfr_set(W.re[AT(i, j)], ar, MPFR_RNDN);
                    mpfr_set(W.im[AT(i, j)], ai, MPFR_RNDN);
                }
            mpfr_clears(ar, ai, (mpfr_ptr)0);
            /* f of a real matrix is real when its imaginary part is rounding */
            bool force_real = false;
            if (real_in) {
                mpfr_t mre, mim;
                mpfr_inits2(bits, mre, mim, (mpfr_ptr)0);
                mpfr_set_zero(mre, 1); mpfr_set_zero(mim, 1);
                for (size_t k = 0; k < nn; k++) {
                    mpfr_abs(tr, W.re[k], MPFR_RNDN); if (mpfr_cmp(tr, mre) > 0) mpfr_set(mre, tr, MPFR_RNDN);
                    mpfr_abs(tr, W.im[k], MPFR_RNDN); if (mpfr_cmp(tr, mim) > 0) mpfr_set(mim, tr, MPFR_RNDN);
                }
                mpfr_mul_2si(mre, mre, 8 - (long)bits, MPFR_RNDN);
                force_real = mpfr_cmp(mim, mre) <= 0;
                mpfr_clears(mre, mim, (mpfr_ptr)0);
            }
            out = mf_mat_to_expr(&W, bits, force_real);
            mf_mat_free(&W);
        }
    }
#undef AT
    mpfr_clears(sr, si, tr, ti, dr, di, den, sep, scale, (mpfr_ptr)0);
    for (size_t k = 0; k < nn; k++) { mpfr_clear(Fr[k]); mpfr_clear(Fi[k]); }
    free(Fr); free(Fi);
    schur_mpfr_core_free(Tr, n); schur_mpfr_core_free(Ti, n);
    schur_mpfr_core_free(Qr, n); schur_mpfr_core_free(Qi, n);
    return out;
}

#endif /* USE_MPFR */
