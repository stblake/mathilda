/*
 * schurdecomp_mpfr.c -- arbitrary-precision SchurDecomposition[m].
 *
 * Handles the standard, real, block-diagonal, un-pivoted case at MPFR
 * precision by reusing the in-house Hessenberg + Francis QR from eigen_direct.c
 * (eigen_schur_real_mpfr), which leaves the real Schur form T in the matrix and
 * the orthogonal Schur vectors Q alongside, with m == Q * T * Q^T.
 *
 * Complex, generalized, RealBlockDiagonalForm -> False, and Pivoting inputs are
 * NOT handled here (a documented gap -- they would need a complex MPFR QR / an
 * MPFR QZ); the dispatcher only calls this kernel for the supported case, and
 * on any failure it falls back to the machine LAPACK kernel.
 *
 * Memory contract: standard builtin ownership (SPEC.md §4).  Never frees the
 * input; every mpfr_t initialised here is cleared on every exit path.
 */

#include "schurdecomp.h"
#include "schurdecomp_internal.h"
#include "eigen.h"
#include "numeric.h"
#include "sym_names.h"
#include "expr.h"
#include "numeric_complex.h"   /* numeric_mpfr_make_complex */

#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

#ifndef USE_MPFR

/* Without MPFR there is no arbitrary-precision Schur kernel: schurdecomp.c only
 * references schur_mpfr_standard_real under #ifdef USE_MPFR, so nothing here is
 * declared or called.  A typedef keeps the translation unit non-empty. */
typedef int schurdecomp_mpfr_no_mpfr_build;

#else

/* Allocate / free a flat array of `count` mpfr_t cells at `bits` precision. */
static mpfr_t* schur_mpfr_alloc(size_t count, mpfr_prec_t bits) {
    mpfr_t* a = (mpfr_t*)malloc(sizeof(mpfr_t) * count);
    for (size_t i = 0; i < count; i++) mpfr_init2(a[i], bits);
    return a;
}

static void schur_mpfr_free(mpfr_t* a, size_t count) {
    for (size_t i = 0; i < count; i++) mpfr_clear(a[i]);
    free(a);
}

/* Load the n x n real matrix (a boxed List-of-Lists -- an NDArray never reaches
 * this arbitrary-precision path) into a row-major mpfr buffer at `bits`.  The
 * imaginary part is discarded (the caller asserts realness via schur_all_real).
 * Returns false, leaking nothing, on a ragged shape or a non-numeric leaf. */
static bool schur_mpfr_load(const Expr* m, int n, mpfr_prec_t bits, mpfr_t* A) {
    if (m->type != EXPR_FUNCTION) return false;
    mpfr_t tmp_im;
    mpfr_init2(tmp_im, bits);
    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (row->type != EXPR_FUNCTION || (int)row->data.function.arg_count != n) {
            mpfr_clear(tmp_im);
            return false;
        }
        for (int j = 0; j < n; j++) {
            bool inexact = false;
            if (!get_approx_mpfr(row->data.function.args[j],
                                 A[(size_t)i * (size_t)n + (size_t)j],
                                 tmp_im, &inexact)) {
                mpfr_clear(tmp_im);
                return false;
            }
        }
    }
    mpfr_clear(tmp_im);
    return true;
}

/* Wrap a row-major n x n mpfr buffer as a nested List of EXPR_MPFR entries. */
static Expr* schur_mpfr_wrap(const mpfr_t* buf, int n) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int j = 0; j < n; j++)
            elems[j] = expr_new_mpfr_copy(buf[(size_t)i * (size_t)n + (size_t)j]);
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), elems, (size_t)n);
        free(elems);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    return out;
}

Expr* schur_mpfr_standard_real(const Expr* m, int n, mpfr_prec_t bits,
                               const SchurOpts* opts) {
    (void)opts;  /* only reached for real / block-diagonal / un-pivoted input */

    size_t nn = (size_t)n * (size_t)n;
    mpfr_t* A = schur_mpfr_alloc(nn, bits);  /* -> Schur form T */
    if (!schur_mpfr_load(m, n, bits, A)) {
        schur_mpfr_free(A, nn);
        return NULL;
    }

    mpfr_t* Q = schur_mpfr_alloc(nn, bits);  /* -> Schur vectors q */
    int status = eigen_schur_real_mpfr(A, n, bits, Q);

    Expr* result = NULL;
    if (status == 0) {
        Expr* q = schur_mpfr_wrap(Q, n);
        Expr* t = schur_mpfr_wrap(A, n);
        Expr* items[2] = { q, t };
        result = expr_new_function(expr_new_symbol(SYM_List), items, 2);
    }

    schur_mpfr_free(A, nn);
    schur_mpfr_free(Q, nn);
    return result;
}

/* ------------------------------------------------------------------------ *
 *  Complex Schur form at arbitrary precision.                               *
 *                                                                           *
 *  m = q . t . ConjugateTranspose[q], t upper triangular, q unitary, every  *
 *  entry at the input precision.  Used for a complex matrix (the real       *
 *  kernel above cannot take one, and the machine fallback dropped a 30-digit *
 *  input to MachinePrecision) and for RealBlockDiagonalForm -> False.       *
 *                                                                           *
 *  Givens rotations throughout: a reduction to upper Hessenberg form, then  *
 *  the single-shift QR iteration with the Wilkinson shift (the eigenvalue   *
 *  of the trailing 2 x 2 nearer its corner), chasing the bulge down the      *
 *  subdiagonal, deflating a subdiagonal entry below eps (|h_kk| + |h_k-1|)  *
 *  and taking an exceptional shift every tenth sweep without progress.      *
 *  The work runs SCHUR_GUARD bits above the input precision and is rounded  *
 *  back on output.                                                          *
 * ------------------------------------------------------------------------ */

#define SCHUR_GUARD 32

typedef struct {
    int n;
    mpfr_t* re;   /* n * n row-major */
    mpfr_t* im;
} SchurCM;

#define SCM(M, i, j) ((size_t)(i) * (size_t)(M)->n + (size_t)(j))

/* |z| */
static void scm_abs(mpfr_t out, mpfr_t re, mpfr_t im) { mpfr_hypot(out, re, im, MPFR_RNDN); }

/* Givens rotation G = [c s; -conj(s) c] with G (a, b)^T = (r, 0)^T. */
static void scm_givens(mpfr_t c, mpfr_t sr, mpfr_t si,
                       mpfr_t ar, mpfr_t ai, mpfr_t br, mpfr_t bi, mpfr_t t1, mpfr_t t2)
{
    scm_abs(t1, ar, ai);                     /* |a| */
    scm_abs(t2, br, bi);                     /* |b| */
    if (mpfr_zero_p(t2)) { mpfr_set_ui(c, 1, MPFR_RNDN); mpfr_set_zero(sr, 1); mpfr_set_zero(si, 1); return; }
    if (mpfr_zero_p(t1)) { mpfr_set_zero(c, 1); mpfr_set_ui(sr, 1, MPFR_RNDN); mpfr_set_zero(si, 1); return; }
    mpfr_t nrm; mpfr_init2(nrm, mpfr_get_prec(c));
    mpfr_hypot(nrm, t1, t2, MPFR_RNDN);       /* sqrt(|a|^2 + |b|^2) */
    mpfr_div(c, t1, nrm, MPFR_RNDN);
    /* s = (a/|a|) conj(b) / nrm */
    mpfr_t pr, pi;
    mpfr_init2(pr, mpfr_get_prec(c)); mpfr_init2(pi, mpfr_get_prec(c));
    /* a conj(b) = (ar br + ai bi) + i (ai br - ar bi) */
    mpfr_mul(pr, ar, br, MPFR_RNDN); mpfr_fma(pr, ai, bi, pr, MPFR_RNDN);
    mpfr_mul(pi, ai, br, MPFR_RNDN); mpfr_fms(pi, ar, bi, pi, MPFR_RNDN); mpfr_neg(pi, pi, MPFR_RNDN);
    mpfr_mul(nrm, nrm, t1, MPFR_RNDN);
    mpfr_div(sr, pr, nrm, MPFR_RNDN);
    mpfr_div(si, pi, nrm, MPFR_RNDN);
    mpfr_clear(pr); mpfr_clear(pi); mpfr_clear(nrm);
}

/* rows i, j <- G (rows i, j), over columns [c0, c1) */
static void scm_rot_rows(SchurCM* M, int i, int j, int c0, int c1,
                         mpfr_t c, mpfr_t sr, mpfr_t si, mpfr_t* w)
{
    /* w[0..5] scratch */
    for (int k = c0; k < c1; k++) {
        size_t x = SCM(M, i, k), y = SCM(M, j, k);
        /* new_x = c x + s y */
        mpfr_mul(w[0], sr, M->re[y], MPFR_RNDN); mpfr_mul(w[1], si, M->im[y], MPFR_RNDN);
        mpfr_sub(w[0], w[0], w[1], MPFR_RNDN);                    /* Re(s y) */
        mpfr_mul(w[2], sr, M->im[y], MPFR_RNDN); mpfr_mul(w[1], si, M->re[y], MPFR_RNDN);
        mpfr_add(w[2], w[2], w[1], MPFR_RNDN);                    /* Im(s y) */
        /* new_y = -conj(s) x + c y ; conj(s) x = (sr xr + si xi) + i (sr xi - si xr) */
        mpfr_mul(w[3], sr, M->re[x], MPFR_RNDN); mpfr_mul(w[1], si, M->im[x], MPFR_RNDN);
        mpfr_add(w[3], w[3], w[1], MPFR_RNDN);
        mpfr_mul(w[4], sr, M->im[x], MPFR_RNDN); mpfr_mul(w[1], si, M->re[x], MPFR_RNDN);
        mpfr_sub(w[4], w[4], w[1], MPFR_RNDN);
        mpfr_mul(w[5], c, M->re[x], MPFR_RNDN); mpfr_add(M->re[x], w[5], w[0], MPFR_RNDN);
        mpfr_mul(w[5], c, M->im[x], MPFR_RNDN); mpfr_add(M->im[x], w[5], w[2], MPFR_RNDN);
        mpfr_mul(w[5], c, M->re[y], MPFR_RNDN); mpfr_sub(M->re[y], w[5], w[3], MPFR_RNDN);
        mpfr_mul(w[5], c, M->im[y], MPFR_RNDN); mpfr_sub(M->im[y], w[5], w[4], MPFR_RNDN);
    }
}

/* columns i, j <- (columns i, j) G^H, over rows [r0, r1) */
static void scm_rot_cols(SchurCM* M, int i, int j, int r0, int r1,
                         mpfr_t c, mpfr_t sr, mpfr_t si, mpfr_t* w)
{
    for (int k = r0; k < r1; k++) {
        size_t x = SCM(M, k, i), y = SCM(M, k, j);
        /* new_x = c x + conj(s) y */
        mpfr_mul(w[0], sr, M->re[y], MPFR_RNDN); mpfr_mul(w[1], si, M->im[y], MPFR_RNDN);
        mpfr_add(w[0], w[0], w[1], MPFR_RNDN);                    /* Re(conj(s) y) */
        mpfr_mul(w[2], sr, M->im[y], MPFR_RNDN); mpfr_mul(w[1], si, M->re[y], MPFR_RNDN);
        mpfr_sub(w[2], w[2], w[1], MPFR_RNDN);                    /* Im(conj(s) y) */
        /* new_y = -s x + c y */
        mpfr_mul(w[3], sr, M->re[x], MPFR_RNDN); mpfr_mul(w[1], si, M->im[x], MPFR_RNDN);
        mpfr_sub(w[3], w[3], w[1], MPFR_RNDN);                    /* Re(s x) */
        mpfr_mul(w[4], sr, M->im[x], MPFR_RNDN); mpfr_mul(w[1], si, M->re[x], MPFR_RNDN);
        mpfr_add(w[4], w[4], w[1], MPFR_RNDN);                    /* Im(s x) */
        mpfr_mul(w[5], c, M->re[x], MPFR_RNDN); mpfr_add(M->re[x], w[5], w[0], MPFR_RNDN);
        mpfr_mul(w[5], c, M->im[x], MPFR_RNDN); mpfr_add(M->im[x], w[5], w[2], MPFR_RNDN);
        mpfr_mul(w[5], c, M->re[y], MPFR_RNDN); mpfr_sub(M->re[y], w[5], w[3], MPFR_RNDN);
        mpfr_mul(w[5], c, M->im[y], MPFR_RNDN); mpfr_sub(M->im[y], w[5], w[4], MPFR_RNDN);
    }
}

/* Apply the rotation built from (a, b) = (H[i][col], H[j][col]) as a unitary
 * similarity on rows / columns i, j of H, and accumulate it into Q. */
static void scm_similarity(SchurCM* H, SchurCM* Q, int i, int j, int col, mpfr_t* w)
{
    int n = H->n;
    mpfr_t c, sr, si;
    mpfr_prec_t p = mpfr_get_prec(w[0]);
    mpfr_init2(c, p); mpfr_init2(sr, p); mpfr_init2(si, p);
    size_t a = SCM(H, i, col), b = SCM(H, j, col);
    scm_givens(c, sr, si, H->re[a], H->im[a], H->re[b], H->im[b], w[6], w[7]);
    scm_rot_rows(H, i, j, 0, n, c, sr, si, w);
    scm_rot_cols(H, i, j, 0, n, c, sr, si, w);
    scm_rot_cols(Q, i, j, 0, n, c, sr, si, w);
    mpfr_set_zero(H->re[b], 1); mpfr_set_zero(H->im[b], 1);
    mpfr_clear(c); mpfr_clear(sr); mpfr_clear(si);
}

/* Shifted QR sweep on the active window [lo, hi] with shift mu: the first
 * rotation is built from (H[lo][lo] - mu, H[lo+1][lo]); the rest chase the
 * bulge H[k+1][k-1] down the subdiagonal. */
static void scm_sweep(SchurCM* H, SchurCM* Q, int lo, int hi,
                      mpfr_t mur, mpfr_t mui, mpfr_t* w)
{
    int n = H->n;
    mpfr_prec_t p = mpfr_get_prec(w[0]);
    mpfr_t c, sr, si, xr, xi;
    mpfr_init2(c, p); mpfr_init2(sr, p); mpfr_init2(si, p);
    mpfr_init2(xr, p); mpfr_init2(xi, p);
    size_t d = SCM(H, lo, lo), b = SCM(H, lo + 1, lo);
    mpfr_sub(xr, H->re[d], mur, MPFR_RNDN);
    mpfr_sub(xi, H->im[d], mui, MPFR_RNDN);
    scm_givens(c, sr, si, xr, xi, H->re[b], H->im[b], w[6], w[7]);
    scm_rot_rows(H, lo, lo + 1, 0, n, c, sr, si, w);
    scm_rot_cols(H, lo, lo + 1, 0, n, c, sr, si, w);
    scm_rot_cols(Q, lo, lo + 1, 0, n, c, sr, si, w);
    for (int k = lo + 1; k < hi; k++) scm_similarity(H, Q, k, k + 1, k - 1, w);
    mpfr_clear(c); mpfr_clear(sr); mpfr_clear(si); mpfr_clear(xr); mpfr_clear(xi);
}

/* The eigenvalue of [[a, b], [c, d]] nearer d: the Wilkinson shift. */
static void scm_wilkinson(SchurCM* H, int hi, mpfr_t mur, mpfr_t mui, mpfr_t* w)
{
    size_t A = SCM(H, hi - 1, hi - 1), B = SCM(H, hi - 1, hi),
           C = SCM(H, hi, hi - 1), D = SCM(H, hi, hi);
    /* h = (a - d)/2;  disc = h^2 + b c;  mu = d + h -+ sqrt(disc), the sign
     * making |h -+ sqrt| smallest i.e. mu nearer d */
    mpfr_t hr, hi_, er, ei, sr, si, t;
    mpfr_prec_t p = mpfr_get_prec(mur);
    mpfr_inits2(p, hr, hi_, er, ei, sr, si, t, (mpfr_ptr)0);
    mpfr_sub(hr, H->re[A], H->re[D], MPFR_RNDN); mpfr_div_2ui(hr, hr, 1, MPFR_RNDN);
    mpfr_sub(hi_, H->im[A], H->im[D], MPFR_RNDN); mpfr_div_2ui(hi_, hi_, 1, MPFR_RNDN);
    /* e = h^2 + b c */
    mpfr_sqr(er, hr, MPFR_RNDN); mpfr_sqr(t, hi_, MPFR_RNDN); mpfr_sub(er, er, t, MPFR_RNDN);
    mpfr_mul(ei, hr, hi_, MPFR_RNDN); mpfr_mul_2ui(ei, ei, 1, MPFR_RNDN);
    mpfr_mul(t, H->re[B], H->re[C], MPFR_RNDN); mpfr_add(er, er, t, MPFR_RNDN);
    mpfr_mul(t, H->im[B], H->im[C], MPFR_RNDN); mpfr_sub(er, er, t, MPFR_RNDN);
    mpfr_mul(t, H->re[B], H->im[C], MPFR_RNDN); mpfr_add(ei, ei, t, MPFR_RNDN);
    mpfr_mul(t, H->im[B], H->re[C], MPFR_RNDN); mpfr_add(ei, ei, t, MPFR_RNDN);
    /* s = sqrt(e): re = sqrt((|e| + er)/2), im = sign(ei) sqrt((|e| - er)/2) */
    mpfr_hypot(t, er, ei, MPFR_RNDN);
    mpfr_add(sr, t, er, MPFR_RNDN); mpfr_div_2ui(sr, sr, 1, MPFR_RNDN);
    if (mpfr_sgn(sr) < 0) mpfr_set_zero(sr, 1);
    mpfr_sqrt(sr, sr, MPFR_RNDN);
    mpfr_sub(si, t, er, MPFR_RNDN); mpfr_div_2ui(si, si, 1, MPFR_RNDN);
    if (mpfr_sgn(si) < 0) mpfr_set_zero(si, 1);
    mpfr_sqrt(si, si, MPFR_RNDN);
    if (mpfr_sgn(ei) < 0) mpfr_neg(si, si, MPFR_RNDN);
    /* pick the sign: compare |h + s| and |h - s|; mu = d + h - (h +- s) ... the
     * eigenvalues are d + h +- s; nearer d is the one with smaller |h +- s| */
    mpfr_t pr, pi, mr, mi;
    mpfr_inits2(p, pr, pi, mr, mi, (mpfr_ptr)0);
    mpfr_add(pr, hr, sr, MPFR_RNDN); mpfr_add(pi, hi_, si, MPFR_RNDN);
    mpfr_sub(mr, hr, sr, MPFR_RNDN); mpfr_sub(mi, hi_, si, MPFR_RNDN);
    mpfr_hypot(t, pr, pi, MPFR_RNDN);
    mpfr_hypot(er, mr, mi, MPFR_RNDN);
    if (mpfr_cmp(t, er) <= 0) { mpfr_add(mur, H->re[D], pr, MPFR_RNDN); mpfr_add(mui, H->im[D], pi, MPFR_RNDN); }
    else                      { mpfr_add(mur, H->re[D], mr, MPFR_RNDN); mpfr_add(mui, H->im[D], mi, MPFR_RNDN); }
    mpfr_clears(hr, hi_, er, ei, sr, si, t, pr, pi, mr, mi, (mpfr_ptr)0);
    (void)w;
}

static Expr* scm_wrap(const SchurCM* M, mpfr_prec_t bits)
{
    int n = M->n;
    mpfr_t r, i;
    mpfr_init2(r, bits); mpfr_init2(i, bits);
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int a = 0; a < n; a++) {
        Expr** el = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int b = 0; b < n; b++) {
            mpfr_set(r, M->re[SCM(M, a, b)], MPFR_RNDN);
            mpfr_set(i, M->im[SCM(M, a, b)], MPFR_RNDN);
            el[b] = mpfr_zero_p(i) ? expr_new_mpfr_copy(r) : numeric_mpfr_make_complex(r, i);
        }
        rows[a] = expr_new_function(expr_new_symbol(SYM_List), el, (size_t)n);
        free(el);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    mpfr_clear(r); mpfr_clear(i);
    return out;
}

static void scm_free(SchurCM* M)
{
    size_t nn = (size_t)M->n * (size_t)M->n;
    schur_mpfr_free(M->re, nn);
    schur_mpfr_free(M->im, nn);
}

bool schur_mpfr_complex_core(const Expr* m, int n, mpfr_prec_t wp,
                             mpfr_t** Tre, mpfr_t** Tim, mpfr_t** Qre, mpfr_t** Qim)
{
    if (m->type != EXPR_FUNCTION) return false;
    size_t nn = (size_t)n * (size_t)n;
    SchurCM H = { n, schur_mpfr_alloc(nn, wp), schur_mpfr_alloc(nn, wp) };
    SchurCM Q = { n, schur_mpfr_alloc(nn, wp), schur_mpfr_alloc(nn, wp) };
    bool ok = true;
    for (int i = 0; i < n && ok; i++) {
        Expr* row = m->data.function.args[i];
        if (row->type != EXPR_FUNCTION || (int)row->data.function.arg_count != n) { ok = false; break; }
        for (int j = 0; j < n; j++) {
            bool inexact = false;
            mpfr_set_zero(H.im[SCM(&H, i, j)], 1);
            if (!get_approx_mpfr(row->data.function.args[j], H.re[SCM(&H, i, j)],
                                 H.im[SCM(&H, i, j)], &inexact)) { ok = false; break; }
            mpfr_set_si(Q.re[SCM(&Q, i, j)], i == j ? 1 : 0, MPFR_RNDN);
            mpfr_set_zero(Q.im[SCM(&Q, i, j)], 1);
        }
    }
    if (!ok) { scm_free(&H); scm_free(&Q); return false; }

    mpfr_t w[8];
    for (int k = 0; k < 8; k++) mpfr_init2(w[k], wp);
    /* Hessenberg: zero H[i][j] for i > j + 1, bottom up, by rotations of
     * rows (i - 1, i) */
    for (int j = 0; j + 2 < n; j++)
        for (int i = n - 1; i >= j + 2; i--) {
            size_t b = SCM(&H, i, j);
            if (mpfr_zero_p(H.re[b]) && mpfr_zero_p(H.im[b])) continue;
            scm_similarity(&H, &Q, i - 1, i, j, w);
        }

    mpfr_t mur, mui, t1, t2, eps;
    mpfr_inits2(wp, mur, mui, t1, t2, eps, (mpfr_ptr)0);
    mpfr_set_ui(eps, 1, MPFR_RNDN);
    mpfr_mul_2si(eps, eps, -(long)wp + 2, MPFR_RNDN);
    int hi = n - 1, iter = 0, total = 0;
    while (hi > 0 && ok) {
        /* deflate: find the active window [lo, hi] */
        int lo = hi;
        while (lo > 0) {
            size_t s = SCM(&H, lo, lo - 1);
            scm_abs(t1, H.re[s], H.im[s]);
            scm_abs(t2, H.re[SCM(&H, lo, lo)], H.im[SCM(&H, lo, lo)]);
            scm_abs(mur, H.re[SCM(&H, lo - 1, lo - 1)], H.im[SCM(&H, lo - 1, lo - 1)]);
            mpfr_add(t2, t2, mur, MPFR_RNDN);
            mpfr_mul(t2, t2, eps, MPFR_RNDN);
            if (mpfr_cmp(t1, t2) <= 0) {
                mpfr_set_zero(H.re[s], 1); mpfr_set_zero(H.im[s], 1);
                break;
            }
            lo--;
        }
        if (lo == hi) { hi--; iter = 0; continue; }
        if (++total > 60 * n) { ok = false; break; }
        if (++iter % 10 == 0) {
            /* exceptional shift: perturb by the subdiagonal magnitude */
            size_t s = SCM(&H, hi, hi - 1);
            scm_abs(t1, H.re[s], H.im[s]);
            mpfr_add(mur, H.re[SCM(&H, hi, hi)], t1, MPFR_RNDN);
            mpfr_set(mui, H.im[SCM(&H, hi, hi)], MPFR_RNDN);
        } else {
            scm_wilkinson(&H, hi, mur, mui, w);
        }
        scm_sweep(&H, &Q, lo, hi, mur, mui, w);
    }
    mpfr_clears(mur, mui, t1, t2, eps, (mpfr_ptr)0);
    for (int k = 0; k < 8; k++) mpfr_clear(w[k]);

    if (!ok) { scm_free(&H); scm_free(&Q); return false; }
    /* t is upper triangular: below the diagonal exactly 0 */
    for (int i = 1; i < n; i++)
        for (int j = 0; j < i; j++) {
            mpfr_set_zero(H.re[SCM(&H, i, j)], 1);
            mpfr_set_zero(H.im[SCM(&H, i, j)], 1);
        }
    *Tre = H.re; *Tim = H.im; *Qre = Q.re; *Qim = Q.im;
    return true;
}

void schur_mpfr_core_free(mpfr_t* a, int n)
{
    schur_mpfr_free(a, (size_t)n * (size_t)n);
}

Expr* schur_mpfr_standard_complex(const Expr* m, int n, mpfr_prec_t bits,
                                  const SchurOpts* opts)
{
    (void)opts;
    SchurCM H = { n, NULL, NULL }, Q = { n, NULL, NULL };
    if (!schur_mpfr_complex_core(m, n, bits + SCHUR_GUARD, &H.re, &H.im, &Q.re, &Q.im))
        return NULL;
    Expr* items[2] = { scm_wrap(&Q, bits), scm_wrap(&H, bits) };
    scm_free(&H); scm_free(&Q);
    return expr_new_function(expr_new_symbol(SYM_List), items, 2);
}

#endif /* USE_MPFR */
