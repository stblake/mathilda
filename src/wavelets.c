/*
 * wavelets.c -- wavelet families and their scaling / wavelet functions.
 *
 *   HaarWavelet[]            the Haar wavelet
 *   DaubechiesWavelet[n]     the Daubechies wavelet with n vanishing moments
 *                            (2n filter taps; DaubechiesWavelet[] = [2])
 *   WaveletPhi[wave, x]      the scaling function phi(x)
 *   WaveletPsi[wave, x]      the wavelet function psi(x)
 *   WaveletFilterCoefficients[wave, "PrimalLowpass" | "PrimalHighpass"]
 *
 * Options of WaveletPhi / WaveletPsi: MaxRecursion -> Automatic (= 8) and
 * WorkingPrecision -> MachinePrecision.
 *
 * Semantics, as in the reference:
 *   - Haar: the exact step functions -- phi = 1 on [0, 1); psi = 1 on
 *     (0, 1/2), -1 on (1/2, 1), 0 elsewhere (psi(0) = psi(1/2) = 0).  A
 *     symbolic x gives the Piecewise.
 *   - Daubechies: phi and psi have no closed form.  The CASCADE gives their
 *     exact values on the dyadic grid k 2^-J (J = MaxRecursion): phi at the
 *     integers is the eigenvector of the two-scale relation phi(x) = 2 Sum_k
 *     h_k phi(2x - k) (normalised Sum phi(k) = 1), the relation itself then
 *     fills each finer level, and psi(x) = 2 Sum_k g_k phi(2x - k) with
 *     g_k = (-1)^k h_{1-k}.  Between grid points the value is the LINEAR
 *     interpolant: WaveletPhi[DaubechiesWavelet[2], 1.3] = 0.34067377391167836...,
 *     the piecewise-linear value at spacing 2^-8, digit for digit the
 *     reference's (verified to 40 digits).  A symbolic x gives
 *     Piecewise[{{InterpolatingFunction[..][x], 0 <= x <= 2n - 1}}, 0] (psi on
 *     [1 - n, n]), the reference's form.
 *   - The filter: h_k with Sum h_k = 1, from the spectral factorisation of the
 *     Daubechies polynomial P(y) = Sum_{k<n} Binomial[n-1+k, k] y^k: each root
 *     y_j of P gives z_j + 1/z_j = 2 - 4 y_j, the root |z_j| < 1 is kept, and
 *     h is the REVERSED coefficient list of (1 + z)^n Prod_j (z - z_j),
 *     normalised -- h = {0.3415, 0.5915, 0.1585, -0.0915} for n = 2, the
 *     reference's ordering.  P's roots by Durand-Kerner in MPFR complex, at the
 *     working precision plus guard bits, so WorkingPrecision -> 50 is 50 digits
 *     throughout.
 *   - Precision: the result has the smaller of WorkingPrecision and the
 *     precision of x (an exact x takes the WorkingPrecision; the default is a
 *     machine number, even for WaveletPhi[w, 1/3]).
 *
 * The cascade grid is memoised per (n, J, precision): NIntegrate or an
 * NDSolve initial condition samples phi thousands of times.
 *
 * Deliberate difference: a List x is mapped element-wise (the reference
 * returns a Piecewise over the whole list).  No NDArray / Compile surface: the
 * first argument is a wavelet object, not a number, and a symbolic x builds an
 * InterpolatingFunction (documented EXEMPT).
 */

#include "wavelets.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "arithmetic.h"
#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "message.h"
#include "numeric.h"
#include "print.h"
#include "sym_names.h"
#include "symtab.h"
#ifdef MATHILDA_THREADS
#include <pthread.h>
#endif

#ifdef USE_MPFR

/* ---- expression helpers (each consumes its arguments) ----------------- */

static Expr* wv_fn(const char* head, Expr** a, size_t n) {
    return expr_new_function(expr_new_symbol(head), a, n);
}
static Expr* wv_fn2(const char* h, Expr* a, Expr* b) { return wv_fn(h, (Expr*[]){ a, b }, 2); }
static Expr* wv_ev(Expr* e) { return eval_and_free(e); }
static size_t wv_len(const Expr* e) { return e->data.function.arg_count; }
static const Expr* wv_at(const Expr* e, size_t i) { return e->data.function.args[i]; }

/* ---- the wavelet object ------------------------------------------------ */

/* 0: not a wavelet; else the number of vanishing moments n (Haar: n = 1 with
 * *haar set). */
static long wv_family(const Expr* w, bool* haar) {
    *haar = false;
    if (head_is(w, SYM_HaarWavelet) && wv_len(w) == 0) { *haar = true; return 1; }
    if (head_is(w, SYM_DaubechiesWavelet)) {
        if (wv_len(w) == 0) return 2;
        const Expr* n = wv_at(w, 0);
        if (wv_len(w) == 1 && n->type == EXPR_INTEGER && n->data.integer >= 1 &&
            n->data.integer <= 512)
            return (long)n->data.integer;
    }
    return 0;
}

/* ---- MPFR complex scratch for the Daubechies roots ---------------------- */

typedef struct { mpfr_t re, im; } wv_c;

static void wv_c_init(wv_c* z, mpfr_prec_t p) { mpfr_inits2(p, z->re, z->im, (mpfr_ptr)0); }
static void wv_c_clear(wv_c* z) { mpfr_clears(z->re, z->im, (mpfr_ptr)0); }
static void wv_c_set(wv_c* z, const wv_c* a) { mpfr_set(z->re, a->re, MPFR_RNDN); mpfr_set(z->im, a->im, MPFR_RNDN); }
/* z = a * b (z may alias a or b); t1, t2 scratch */
static void wv_c_mul(wv_c* z, const wv_c* a, const wv_c* b, mpfr_t t1, mpfr_t t2) {
    mpfr_mul(t1, a->re, b->re, MPFR_RNDN);
    mpfr_mul(t2, a->im, b->im, MPFR_RNDN);
    mpfr_sub(t1, t1, t2, MPFR_RNDN);
    mpfr_mul(t2, a->re, b->im, MPFR_RNDN);
    mpfr_fma(t2, a->im, b->re, t2, MPFR_RNDN);
    mpfr_set(z->re, t1, MPFR_RNDN);
    mpfr_set(z->im, t2, MPFR_RNDN);
}
/* z = a / b */
static void wv_c_div(wv_c* z, const wv_c* a, const wv_c* b, mpfr_t t1, mpfr_t t2, mpfr_t t3) {
    mpfr_sqr(t3, b->re, MPFR_RNDN);
    mpfr_fma(t3, b->im, b->im, t3, MPFR_RNDN);
    mpfr_mul(t1, a->re, b->re, MPFR_RNDN);
    mpfr_fma(t1, a->im, b->im, t1, MPFR_RNDN);
    mpfr_mul(t2, a->im, b->re, MPFR_RNDN);
    mpfr_mul(z->re, a->re, b->im, MPFR_RNDN);
    mpfr_sub(t2, t2, z->re, MPFR_RNDN);
    mpfr_div(z->re, t1, t3, MPFR_RNDN);
    mpfr_div(z->im, t2, t3, MPFR_RNDN);
}

/* The Daubechies lowpass filter h[0..2n-1] (Sum h = 1) at precision p. */
static bool wv_daubechies_filter(long n, mpfr_prec_t p, mpfr_t* h) {
    long L = 2 * n;
    mpfr_prec_t wp = p + 64 + (mpfr_prec_t)(4 * n);
    mpfr_t t1, t2, t3;
    mpfr_inits2(wp, t1, t2, t3, (mpfr_ptr)0);
    long m = n - 1;                                   /* deg P */
    wv_c* zr = NULL;                                  /* the kept z-roots */
    bool ok = true;
    if (m > 0) {
        /* P(y) = Sum_{k<n} C(n-1+k, k) y^k, made monic for Durand-Kerner. */
        mpfr_t* a = malloc(sizeof(mpfr_t) * (size_t)(m + 1));
        for (long k = 0; k <= m; k++) {
            mpfr_init2(a[k], wp);
            mpz_t bz;
            mpz_init(bz);
            mpz_bin_uiui(bz, (unsigned long)(n - 1 + k), (unsigned long)k);
            mpfr_set_z(a[k], bz, MPFR_RNDN);
            mpz_clear(bz);
        }
        for (long k = 0; k <= m; k++) mpfr_div(a[k], a[k], a[m], MPFR_RNDN);
        wv_c* y = malloc(sizeof(wv_c) * (size_t)m);
        wv_c num, den, d, tz;
        wv_c_init(&num, wp); wv_c_init(&den, wp); wv_c_init(&d, wp); wv_c_init(&tz, wp);
        /* starting points on a circle of the Cauchy bound's radius, rotated
         * off the axes so no two coincide or sit on the real line */
        mpfr_t R; mpfr_init2(R, wp);
        mpfr_set_ui(R, 0, MPFR_RNDN);
        for (long k = 0; k < m; k++) {
            mpfr_abs(t1, a[k], MPFR_RNDN);
            if (mpfr_cmp(t1, R) > 0) mpfr_set(R, t1, MPFR_RNDN);
        }
        mpfr_add_ui(R, R, 1, MPFR_RNDN);
        mpfr_mul_d(R, R, 0.5, MPFR_RNDN);
        for (long j = 0; j < m; j++) {
            wv_c_init(&y[j], wp);
            double ang = 2.0 * 3.14159265358979323846 * ((double)j + 0.25) / (double)m + 0.4;
            mpfr_mul_d(y[j].re, R, cos(ang), MPFR_RNDN);
            mpfr_mul_d(y[j].im, R, sin(ang), MPFR_RNDN);
        }
        bool conv = false;
        for (int it = 0; it < 4000 && !conv; it++) {
            conv = true;
            for (long j = 0; j < m; j++) {
                /* num = P(y_j) (Horner), den = Prod_{i != j} (y_j - y_i) */
                mpfr_set_ui(num.re, 1, MPFR_RNDN); mpfr_set_ui(num.im, 0, MPFR_RNDN);
                for (long k = m - 1; k >= 0; k--) {
                    wv_c_mul(&num, &num, &y[j], t1, t2);
                    mpfr_add(num.re, num.re, a[k], MPFR_RNDN);
                }
                mpfr_set_ui(den.re, 1, MPFR_RNDN); mpfr_set_ui(den.im, 0, MPFR_RNDN);
                for (long i = 0; i < m; i++) {
                    if (i == j) continue;
                    mpfr_sub(tz.re, y[j].re, y[i].re, MPFR_RNDN);
                    mpfr_sub(tz.im, y[j].im, y[i].im, MPFR_RNDN);
                    wv_c_mul(&den, &den, &tz, t1, t2);
                }
                wv_c_div(&d, &num, &den, t1, t2, t3);
                mpfr_sub(y[j].re, y[j].re, d.re, MPFR_RNDN);
                mpfr_sub(y[j].im, y[j].im, d.im, MPFR_RNDN);
                /* converged when the correction is below 2^-(wp-8) |y_j| */
                mpfr_hypot(t1, d.re, d.im, MPFR_RNDN);
                mpfr_hypot(t2, y[j].re, y[j].im, MPFR_RNDN);
                mpfr_mul_2si(t2, t2, -(long)(wp - 8), MPFR_RNDN);
                if (mpfr_cmp(t1, t2) > 0) conv = false;
            }
        }
        ok = conv;
        /* z + 1/z = 2 - 4y: z = b +- Sqrt[b^2 - 1], b = 1 - 2y; keep |z| < 1 */
        zr = malloc(sizeof(wv_c) * (size_t)m);
        wv_c b, s;
        wv_c_init(&b, wp); wv_c_init(&s, wp);
        for (long j = 0; j < m; j++) {
            wv_c_init(&zr[j], wp);
            mpfr_mul_si(b.re, y[j].re, -2, MPFR_RNDN); mpfr_add_ui(b.re, b.re, 1, MPFR_RNDN);
            mpfr_mul_si(b.im, y[j].im, -2, MPFR_RNDN);
            wv_c_mul(&s, &b, &b, t1, t2);
            mpfr_sub_ui(s.re, s.re, 1, MPFR_RNDN);
            /* complex sqrt of s */
            mpfr_hypot(t1, s.re, s.im, MPFR_RNDN);
            mpfr_add(t2, t1, s.re, MPFR_RNDN); mpfr_div_2ui(t2, t2, 1, MPFR_RNDN);
            if (mpfr_sgn(t2) < 0) mpfr_set_ui(t2, 0, MPFR_RNDN);
            mpfr_sqrt(t2, t2, MPFR_RNDN);                     /* Re sqrt */
            mpfr_sub(t3, t1, s.re, MPFR_RNDN); mpfr_div_2ui(t3, t3, 1, MPFR_RNDN);
            if (mpfr_sgn(t3) < 0) mpfr_set_ui(t3, 0, MPFR_RNDN);
            mpfr_sqrt(t3, t3, MPFR_RNDN);                     /* |Im sqrt| */
            if (mpfr_sgn(s.im) < 0) mpfr_neg(t3, t3, MPFR_RNDN);
            /* z1 = b + sq; z2 = b - sq: keep the one inside the unit circle */
            mpfr_add(zr[j].re, b.re, t2, MPFR_RNDN);
            mpfr_add(zr[j].im, b.im, t3, MPFR_RNDN);
            mpfr_hypot(t1, zr[j].re, zr[j].im, MPFR_RNDN);
            if (mpfr_cmp_ui(t1, 1) > 0) {
                mpfr_sub(zr[j].re, b.re, t2, MPFR_RNDN);
                mpfr_sub(zr[j].im, b.im, t3, MPFR_RNDN);
            }
        }
        wv_c_clear(&b); wv_c_clear(&s);
        for (long j = 0; j < m; j++) wv_c_clear(&y[j]);
        free(y);
        wv_c_clear(&num); wv_c_clear(&den); wv_c_clear(&d); wv_c_clear(&tz);
        mpfr_clear(R);
        for (long k = 0; k <= m; k++) mpfr_clear(a[k]);
        free(a);
    }
    /* q(z) = (1 + z)^n Prod_j (z - z_j), complex coefficients c[0..L-1] */
    wv_c* c = malloc(sizeof(wv_c) * (size_t)L);
    for (long k = 0; k < L; k++) { wv_c_init(&c[k], wp); mpfr_set_ui(c[k].re, 0, MPFR_RNDN); mpfr_set_ui(c[k].im, 0, MPFR_RNDN); }
    mpfr_set_ui(c[0].re, 1, MPFR_RNDN);
    long deg = 0;
    wv_c tmp; wv_c_init(&tmp, wp);
    for (long f = 0; f < n + m; f++) {
        /* multiply by (z - r), r = -1 for the first n factors, z_j after */
        wv_c r;
        wv_c_init(&r, wp);
        if (f < n) { mpfr_set_si(r.re, -1, MPFR_RNDN); mpfr_set_ui(r.im, 0, MPFR_RNDN); }
        else wv_c_set(&r, &zr[f - n]);
        for (long k = deg + 1; k >= 0; k--) {
            /* c[k] = c[k-1] - r c[k] */
            wv_c_mul(&tmp, &r, &c[k], t1, t2);
            if (k > 0) {
                mpfr_sub(c[k].re, c[k - 1].re, tmp.re, MPFR_RNDN);
                mpfr_sub(c[k].im, c[k - 1].im, tmp.im, MPFR_RNDN);
            } else {
                mpfr_neg(c[k].re, tmp.re, MPFR_RNDN);
                mpfr_neg(c[k].im, tmp.im, MPFR_RNDN);
            }
        }
        deg++;
        wv_c_clear(&r);
    }
    wv_c_clear(&tmp);
    /* h_k = Re c_{L-1-k} / Sum Re c (the imaginary parts cancel in conjugate
     * pairs) */
    mpfr_set_ui(t3, 0, MPFR_RNDN);
    for (long k = 0; k < L; k++) mpfr_add(t3, t3, c[k].re, MPFR_RNDN);
    if (mpfr_zero_p(t3)) ok = false;
    for (long k = 0; k < L; k++) {
        mpfr_init2(h[k], p);
        mpfr_div(h[k], c[L - 1 - k].re, t3, MPFR_RNDN);
    }
    for (long k = 0; k < L; k++) wv_c_clear(&c[k]);
    free(c);
    if (zr) { for (long j = 0; j < m; j++) wv_c_clear(&zr[j]); free(zr); }
    mpfr_clears(t1, t2, t3, (mpfr_ptr)0);
    return ok;
}

/* ---- the cascade grid --------------------------------------------------- */

typedef struct {
    long n, J;
    mpfr_prec_t prec;
    long nphi, npsi;      /* grid sizes: (2n-1) 2^J + 1 each */
    mpfr_t* phi;          /* phi(k 2^-J), k = 0..nphi-1, on [0, 2n-1] */
    mpfr_t* psi;          /* psi((1-n) + k 2^-J), on [1-n, n] */
    double* dphi;         /* the same rounded to double */
    double* dpsi;
} WvGrid;

static void wv_grid_free(WvGrid* g) {
    if (!g) return;
    for (long k = 0; k < g->nphi; k++) mpfr_clear(g->phi[k]);
    for (long k = 0; k < g->npsi; k++) mpfr_clear(g->psi[k]);
    free(g->phi); free(g->psi); free(g->dphi); free(g->dpsi);
    free(g);
}

/* Solve the (L-1) x (L-1) system for phi at the interior integers. */
static bool wv_integer_values(long n, mpfr_t* h, mpfr_prec_t p, mpfr_t* phi_int) {
    long L = 2 * n - 1;                     /* support [0, L] */
    if (n == 1) {                           /* Haar cascade: phi(0) = 1, phi(1) = 0 */
        mpfr_set_ui(phi_int[0], 1, MPFR_RNDN);
        mpfr_set_ui(phi_int[1], 0, MPFR_RNDN);
        return true;
    }
    long N = L - 1;                         /* unknowns phi(1..L-1) */
    mpfr_t* A = malloc(sizeof(mpfr_t) * (size_t)(N * (N + 1)));
    for (long i = 0; i < N * (N + 1); i++) { mpfr_init2(A[i], p); mpfr_set_ui(A[i], 0, MPFR_RNDN); }
    /* rows 0..N-2: phi(i) - 2 Sum_j h_{2i-j} phi(j) = 0; last row: Sum phi = 1 */
    for (long r = 0; r < N - 1; r++) {
        long i = r + 1;
        for (long cidx = 0; cidx < N; cidx++) {
            long j = cidx + 1, k = 2 * i - j;
            if (k >= 0 && k < 2 * n) mpfr_mul_2ui(A[r * (N + 1) + cidx], h[k], 1, MPFR_RNDN);
            mpfr_neg(A[r * (N + 1) + cidx], A[r * (N + 1) + cidx], MPFR_RNDN);
            if (j == i) mpfr_add_ui(A[r * (N + 1) + cidx], A[r * (N + 1) + cidx], 1, MPFR_RNDN);
        }
    }
    for (long cidx = 0; cidx < N; cidx++) mpfr_set_ui(A[(N - 1) * (N + 1) + cidx], 1, MPFR_RNDN);
    mpfr_set_ui(A[(N - 1) * (N + 1) + N], 1, MPFR_RNDN);
    /* Gaussian elimination with partial pivoting */
    mpfr_t t, f;
    mpfr_inits2(p, t, f, (mpfr_ptr)0);
    bool ok = true;
    for (long col = 0; col < N && ok; col++) {
        long piv = col;
        for (long r = col + 1; r < N; r++) {
            mpfr_abs(t, A[r * (N + 1) + col], MPFR_RNDN);
            mpfr_abs(f, A[piv * (N + 1) + col], MPFR_RNDN);
            if (mpfr_cmp(t, f) > 0) piv = r;
        }
        if (mpfr_zero_p(A[piv * (N + 1) + col])) { ok = false; break; }
        if (piv != col)
            for (long cc = 0; cc <= N; cc++) mpfr_swap(A[piv * (N + 1) + cc], A[col * (N + 1) + cc]);
        for (long r = 0; r < N; r++) {
            if (r == col || mpfr_zero_p(A[r * (N + 1) + col])) continue;
            mpfr_div(f, A[r * (N + 1) + col], A[col * (N + 1) + col], MPFR_RNDN);
            for (long cc = col; cc <= N; cc++) {
                mpfr_mul(t, f, A[col * (N + 1) + cc], MPFR_RNDN);
                mpfr_sub(A[r * (N + 1) + cc], A[r * (N + 1) + cc], t, MPFR_RNDN);
            }
        }
    }
    mpfr_set_ui(phi_int[0], 0, MPFR_RNDN);
    mpfr_set_ui(phi_int[L], 0, MPFR_RNDN);
    for (long r = 0; r < N && ok; r++)
        mpfr_div(phi_int[r + 1], A[r * (N + 1) + N], A[r * (N + 1) + r], MPFR_RNDN);
    mpfr_clears(t, f, (mpfr_ptr)0);
    for (long i = 0; i < N * (N + 1); i++) mpfr_clear(A[i]);
    free(A);
    return ok;
}

static WvGrid* wv_grid_build(long n, long J, mpfr_prec_t prec) {
    mpfr_prec_t p = prec + 32;
    long taps = 2 * n, L = 2 * n - 1;
    mpfr_t* h = malloc(sizeof(mpfr_t) * (size_t)taps);
    bool ok;
    if (n == 1) {
        for (long k = 0; k < 2; k++) { mpfr_init2(h[k], p); mpfr_set_d(h[k], 0.5, MPFR_RNDN); }
        ok = true;
    } else {
        ok = wv_daubechies_filter(n, p, h);
    }
    if (!ok) { for (long k = 0; k < taps; k++) mpfr_clear(h[k]); free(h); return NULL; }

    WvGrid* g = calloc(1, sizeof(WvGrid));
    g->n = n; g->J = J; g->prec = prec;
    long S = (long)1 << J;
    g->nphi = L * S + 1;
    g->npsi = L * S + 1;
    g->phi = malloc(sizeof(mpfr_t) * (size_t)g->nphi);
    for (long k = 0; k < g->nphi; k++) { mpfr_init2(g->phi[k], p); mpfr_set_ui(g->phi[k], 0, MPFR_RNDN); }
    mpfr_t* pint = malloc(sizeof(mpfr_t) * (size_t)(L + 1));
    for (long k = 0; k <= L; k++) mpfr_init2(pint[k], p);
    ok = wv_integer_values(n, h, p, pint);
    for (long k = 0; k <= L; k++) { if (ok) mpfr_set(g->phi[k * S], pint[k], MPFR_RNDN); mpfr_clear(pint[k]); }
    free(pint);
    mpfr_t acc, t;
    mpfr_inits2(p, acc, t, (mpfr_ptr)0);
    /* refine: level j fills the odd multiples of 2^(J-j) */
    for (long j = 1; j <= J && ok; j++) {
        long step = (long)1 << (J - j);
        for (long m = step; m < g->nphi; m += 2 * step) {
            mpfr_set_ui(acc, 0, MPFR_RNDN);
            for (long k = 0; k < taps; k++) {
                long idx = 2 * m - k * S;
                if (idx < 0 || idx >= g->nphi) continue;
                mpfr_mul(t, h[k], g->phi[idx], MPFR_RNDN);
                mpfr_add(acc, acc, t, MPFR_RNDN);
            }
            mpfr_mul_2ui(g->phi[m], acc, 1, MPFR_RNDN);
        }
    }
    /* psi((1-n) + m 2^-J) = 2 Sum_{k=2-2n}^{1} (-1)^k h_{1-k} phi(2x - k) */
    g->psi = malloc(sizeof(mpfr_t) * (size_t)g->npsi);
    for (long m = 0; m < g->npsi; m++) {
        mpfr_init2(g->psi[m], p);
        mpfr_set_ui(acc, 0, MPFR_RNDN);
        long xi = (1 - n) * S + m;                       /* x in units of 2^-J */
        for (long k = 2 - 2 * n; k <= 1 && ok; k++) {
            long idx = 2 * xi - k * S;                   /* (2x - k) in units of 2^-J */
            if (idx < 0 || idx >= g->nphi) continue;
            mpfr_mul(t, h[1 - k], g->phi[idx], MPFR_RNDN);
            if (k & 1) mpfr_sub(acc, acc, t, MPFR_RNDN);
            else       mpfr_add(acc, acc, t, MPFR_RNDN);
        }
        mpfr_mul_2ui(g->psi[m], acc, 1, MPFR_RNDN);
    }
    mpfr_clears(acc, t, (mpfr_ptr)0);
    for (long k = 0; k < taps; k++) mpfr_clear(h[k]);
    free(h);
    if (!ok) { wv_grid_free(g); return NULL; }
    g->dphi = malloc(sizeof(double) * (size_t)g->nphi);
    g->dpsi = malloc(sizeof(double) * (size_t)g->npsi);
    for (long k = 0; k < g->nphi; k++) g->dphi[k] = mpfr_get_d(g->phi[k], MPFR_RNDN);
    for (long k = 0; k < g->npsi; k++) g->dpsi[k] = mpfr_get_d(g->psi[k], MPFR_RNDN);
    return g;
}

/* One-entry memo of the last grid (the common use samples one wavelet). */
#ifdef MATHILDA_THREADS
static pthread_mutex_t g_wv_lock = PTHREAD_MUTEX_INITIALIZER;
#define WV_LOCK()   pthread_mutex_lock(&g_wv_lock)
#define WV_UNLOCK() pthread_mutex_unlock(&g_wv_lock)
#else
#define WV_LOCK()   ((void)0)
#define WV_UNLOCK() ((void)0)
#endif
static WvGrid* g_wv_memo;

/* The grid for (n, J, prec); held under the lock by the caller until
 * wv_grid_release. */
static WvGrid* wv_grid_acquire(long n, long J, mpfr_prec_t prec) {
    WV_LOCK();
    if (g_wv_memo && g_wv_memo->n == n && g_wv_memo->J == J && g_wv_memo->prec == prec)
        return g_wv_memo;
    WvGrid* g = wv_grid_build(n, J, prec);
    if (g) { wv_grid_free(g_wv_memo); g_wv_memo = g; }
    if (!g) WV_UNLOCK();
    return g;
}
static void wv_grid_release(void) { WV_UNLOCK(); }

/* ---- evaluation ---------------------------------------------------------- */

/* Interpolated value at x (machine): grid value array, its origin x0 (in
 * units of 1), spacing 2^-J. */
static double wv_interp_d(const double* v, long nv, long x0, long J, double x) {
    double t = ldexp(x - (double)x0, (int)J);
    if (!(t >= 0.0) || t > (double)(nv - 1)) return 0.0;
    long i = (long)floor(t);
    if (i >= nv - 1) return v[nv - 1];
    double f = t - (double)i;
    return v[i] * (1.0 - f) + v[i + 1] * f;
}

static Expr* wv_interp_mpfr(mpfr_t* v, long nv, long x0, long J, const mpfr_t x, mpfr_prec_t outp) {
    mpfr_t t, f, a, b;
    mpfr_inits2(mpfr_get_prec(v[0]) + 32, t, f, a, b, (mpfr_ptr)0);
    mpfr_sub_si(t, x, x0, MPFR_RNDN);
    mpfr_mul_2si(t, t, J, MPFR_RNDN);
    Expr* out = expr_new_mpfr_bits(outp);
    if (mpfr_sgn(t) >= 0 && mpfr_cmp_si(t, nv - 1) <= 0) {
        mpfr_floor(f, t);
        long i = mpfr_get_si(f, MPFR_RNDN);
        if (i >= nv - 1) {
            mpfr_set(out->data.mpfr, v[nv - 1], MPFR_RNDN);
        } else {
            mpfr_sub(f, t, f, MPFR_RNDN);                   /* fraction */
            mpfr_ui_sub(a, 1, f, MPFR_RNDN);
            mpfr_mul(a, a, v[i], MPFR_RNDN);
            mpfr_mul(b, f, v[i + 1], MPFR_RNDN);
            mpfr_add(out->data.mpfr, a, b, MPFR_RNDN);
        }
    }
    mpfr_clears(t, f, a, b, (mpfr_ptr)0);
    return out;
}

/* Haar at a numeric x (as a double, or the comparison sign at MPFR): the
 * exact step values. */
static int wv_haar_value(bool psi, int cmp0, int cmp_half, int cmp1) {
    if (!psi) return (cmp0 >= 0 && cmp1 < 0) ? 1 : 0;
    if (cmp0 > 0 && cmp_half < 0) return 1;
    if (cmp_half > 0 && cmp1 < 0) return -1;
    return 0;
}

/* Symbolic x: the Piecewise forms. */
static Expr* wv_symbolic(bool psi, bool haar, long n, long J, mpfr_prec_t prec, bool machine,
                         const Expr* x) {
    if (haar) {
        Expr* half = make_rational(1, 2);
        if (!psi) {
            Expr* cond = wv_fn("Inequality", (Expr*[]){ expr_new_integer(0), expr_new_symbol(SYM_LessEqual),
                                   expr_copy((Expr*)x), expr_new_symbol(SYM_Less), expr_new_integer(1) }, 5);
            expr_free(half);
            Expr* pair = wv_fn2(SYM_List, expr_new_integer(1), cond);
            return wv_ev(wv_fn2(SYM_Piecewise, wv_fn(SYM_List, (Expr*[]){ pair }, 1), expr_new_integer(0)));
        }
        Expr* c1 = wv_fn(SYM_Less, (Expr*[]){ expr_copy(half), expr_copy((Expr*)x), expr_new_integer(1) }, 3);
        Expr* c2 = wv_fn(SYM_Less, (Expr*[]){ expr_new_integer(0), expr_copy((Expr*)x), half }, 3);
        Expr* p1 = wv_fn2(SYM_List, expr_new_integer(-1), c1);
        Expr* p2 = wv_fn2(SYM_List, expr_new_integer(1), c2);
        return wv_ev(wv_fn2(SYM_Piecewise, wv_fn2(SYM_List, p1, p2), expr_new_integer(0)));
    }
    WvGrid* g = wv_grid_acquire(n, J, prec);
    if (!g) return NULL;
    long nv = psi ? g->npsi : g->nphi;
    long x0 = psi ? 1 - n : 0;
    Expr** rows = malloc(sizeof(Expr*) * (size_t)nv);
    for (long k = 0; k < nv; k++) {
        Expr* xv;
        Expr* yv;
        if (machine) {
            xv = expr_new_real((double)x0 + ldexp((double)k, (int)-J));
            yv = expr_new_real(psi ? g->dpsi[k] : g->dphi[k]);
        } else {
            xv = expr_new_mpfr_bits(prec);
            mpfr_set_si(xv->data.mpfr, x0, MPFR_RNDN);
            mpfr_t d; mpfr_init2(d, prec);
            mpfr_set_si(d, k, MPFR_RNDN);
            mpfr_mul_2si(d, d, -J, MPFR_RNDN);
            mpfr_add(xv->data.mpfr, xv->data.mpfr, d, MPFR_RNDN);
            mpfr_clear(d);
            yv = expr_new_mpfr_bits(prec);
            mpfr_set(yv->data.mpfr, psi ? g->psi[k] : g->phi[k], MPFR_RNDN);
        }
        rows[k] = wv_fn2(SYM_List, xv, yv);
    }
    wv_grid_release();
    Expr* data = wv_fn(SYM_List, rows, (size_t)nv);
    free(rows);
    Expr* ifun = wv_ev(wv_fn2(SYM_Interpolation, data,
                    wv_fn2(SYM_Rule, expr_new_symbol(SYM_InterpolationOrder), expr_new_integer(1))));
    Expr* app = expr_new_function(ifun, (Expr*[]){ expr_copy((Expr*)x) }, 1);
    long x1 = psi ? n : 2 * n - 1;
    Expr* cond = wv_fn(SYM_LessEqual, (Expr*[]){ expr_new_integer(x0), expr_copy((Expr*)x),
                                                expr_new_integer(x1) }, 3);
    Expr* pair = wv_fn2(SYM_List, app, cond);
    return wv_fn2(SYM_Piecewise, wv_fn(SYM_List, (Expr*[]){ pair }, 1), expr_new_integer(0));
}

/* The value at one x.  prec: the result's mantissa width (machine when
 * `machine`). */
static Expr* wv_point(bool psi, bool haar, long n, long J, mpfr_prec_t wprec, bool wp_machine,
                      const Expr* x) {
    bool x_machine = x->type == EXPR_REAL;
    bool x_mpfr = x->type == EXPR_MPFR;
    if (!x_machine && !x_mpfr && !expr_is_numeric_quantity(x))
        return wv_symbolic(psi, haar, n, J, wprec, wp_machine, x);
    Expr *re, *im;
    if (is_complex((Expr*)x, &re, &im)) return NULL;
    /* the result precision: min(WorkingPrecision, precision of x) */
    bool machine = wp_machine || x_machine;
    mpfr_prec_t prec = wprec;
    if (x_mpfr && mpfr_get_prec(x->data.mpfr) < prec) prec = mpfr_get_prec(x->data.mpfr);
    mpfr_t xv;
    mpfr_init2(xv, (wp_machine ? 64 : wprec) + 32);
    if (x_machine) mpfr_set_d(xv, x->data.real, MPFR_RNDN);
    else if (x_mpfr) mpfr_set(xv, x->data.mpfr, MPFR_RNDN);
    else {
        NumericSpec spec;
        spec.mode = NUMERIC_MODE_MPFR;
        spec.bits = (long)mpfr_get_prec(xv);
        Expr* v = numericalize((Expr*)x, spec);
        bool ok = v && (v->type == EXPR_MPFR || v->type == EXPR_REAL);
        if (ok) {
            if (v->type == EXPR_MPFR) mpfr_set(xv, v->data.mpfr, MPFR_RNDN);
            else mpfr_set_d(xv, v->data.real, MPFR_RNDN);
        }
        expr_free(v);
        if (!ok) { mpfr_clear(xv); return NULL; }
    }
    Expr* out = NULL;
    if (haar) {
        int c0 = mpfr_sgn(xv), c1 = mpfr_cmp_ui(xv, 1);
        int ch = mpfr_cmp_d(xv, 0.5);
        int v = wv_haar_value(psi, c0, ch, c1);
        if (machine) out = expr_new_real((double)v);
        else { out = expr_new_mpfr_bits(prec); mpfr_set_si(out->data.mpfr, v, MPFR_RNDN); }
    } else {
        WvGrid* g = wv_grid_acquire(n, J, machine ? 53 : wprec);
        if (g) {
            long nv = psi ? g->npsi : g->nphi;
            long x0 = psi ? 1 - n : 0;
            if (machine) {
                out = expr_new_real(wv_interp_d(psi ? g->dpsi : g->dphi, nv, x0, J,
                                                mpfr_get_d(xv, MPFR_RNDN)));
            } else {
                out = wv_interp_mpfr(psi ? g->psi : g->phi, nv, x0, J, xv, prec);
            }
            wv_grid_release();
        }
    }
    mpfr_clear(xv);
    return out;
}

static Expr* wv_builtin(Expr* res, bool psi) {
    const char* head = psi ? "WaveletPsi" : "WaveletPhi";
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = wv_len(res);
    if (argc < 2) return argc == 1 ? NULL : builtin_arg_error(head, argc, 1, 2);
    bool haar;
    long n = wv_family(wv_at(res, 0), &haar);
    if (n == 0) return NULL;
    /* options */
    long J = 8;
    bool wp_machine = true;
    mpfr_prec_t wprec = 53;
    for (size_t i = 2; i < argc; i++) {
        const Expr* o = wv_at(res, i);
        if (!(head_is(o, SYM_Rule) || head_is(o, SYM_RuleDelayed)) || wv_len(o) != 2) return NULL;
        const Expr* k = wv_at(o, 0);
        const Expr* v = wv_at(o, 1);
        if (k->type != EXPR_SYMBOL) return NULL;
        if (k->data.symbol.name == SYM_MaxRecursion) {
            if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_Automatic) continue;
            if (v->type != EXPR_INTEGER || v->data.integer < 0 || v->data.integer > 24) return NULL;
            J = (long)v->data.integer;
        } else if (k->data.symbol.name == SYM_WorkingPrecision) {
            if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_MachinePrecision) continue;
            double d;
            if (v->type == EXPR_INTEGER) d = (double)v->data.integer;
            else if (v->type == EXPR_REAL) d = v->data.real;
            else return NULL;
            if (!(d > 0) || d > 100000) return NULL;
            wp_machine = false;
            wprec = (mpfr_prec_t)numeric_digits_to_bits(d);
        } else {
            return NULL;
        }
    }
    const Expr* x = wv_at(res, 1);
    if (head_is(x, SYM_List)) {
        size_t m = wv_len(x);
        Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
        for (size_t i = 0; i < m; i++) {
            out[i] = wv_point(psi, haar, n, J, wprec, wp_machine, wv_at(x, i));
            if (!out[i]) {
                for (size_t j = 0; j < i; j++) expr_free(out[j]);
                free(out);
                return NULL;
            }
        }
        Expr* r = wv_fn(SYM_List, out, m);
        free(out);
        return r;
    }
    return wv_point(psi, haar, n, J, wprec, wp_machine, x);
}

Expr* builtin_waveletphi(Expr* res) { return wv_builtin(res, false); }
Expr* builtin_waveletpsi(Expr* res) { return wv_builtin(res, true); }

/* WaveletFilterCoefficients[wave, "PrimalLowpass" | "PrimalHighpass"]:
 * {{k, h_k}, ...}, machine (or WorkingPrecision -> p). */
Expr* builtin_waveletfiltercoefficients(Expr* res) {
    if (res->type != EXPR_FUNCTION || wv_len(res) < 1) return NULL;
    bool haar;
    long n = wv_family(wv_at(res, 0), &haar);
    if (n == 0) return NULL;
    bool high = false;
    bool wp_machine = true;
    mpfr_prec_t wprec = 53;
    for (size_t i = 1; i < wv_len(res); i++) {
        const Expr* a = wv_at(res, i);
        if (a->type == EXPR_STRING) {
            if (strcmp(a->data.string, "PrimalLowpass") == 0 || strcmp(a->data.string, "DualLowpass") == 0) high = false;
            else if (strcmp(a->data.string, "PrimalHighpass") == 0 || strcmp(a->data.string, "DualHighpass") == 0) high = true;
            else return NULL;
            continue;
        }
        if (head_is(a, SYM_Rule) && wv_len(a) == 2 && wv_at(a, 0)->type == EXPR_SYMBOL &&
            wv_at(a, 0)->data.symbol.name == SYM_WorkingPrecision) {
            const Expr* v = wv_at(a, 1);
            if (v->type == EXPR_SYMBOL && v->data.symbol.name == SYM_MachinePrecision) continue;
            double d;
            if (v->type == EXPR_INTEGER) d = (double)v->data.integer;
            else if (v->type == EXPR_REAL) d = v->data.real;
            else return NULL;
            if (!(d > 0) || d > 100000) return NULL;
            wp_machine = false;
            wprec = (mpfr_prec_t)numeric_digits_to_bits(d);
            continue;
        }
        return NULL;
    }
    long taps = 2 * n;
    mpfr_prec_t p = (wp_machine ? 53 : wprec) + 32;
    mpfr_t* h = malloc(sizeof(mpfr_t) * (size_t)taps);
    bool ok = true;
    if (n == 1) { for (long k = 0; k < 2; k++) { mpfr_init2(h[k], p); mpfr_set_d(h[k], 0.5, MPFR_RNDN); } }
    else ok = wv_daubechies_filter(n, p, h);
    Expr* r = NULL;
    if (ok) {
        Expr** rows = malloc(sizeof(Expr*) * (size_t)taps);
        for (long i = 0; i < taps; i++) {
            /* lowpass: k = i, h_k; highpass: k = 2 - 2n + i, (-1)^k h_{1-k} */
            long k = high ? 2 - 2 * n + i : i;
            mpfr_srcptr src = high ? h[1 - k] : h[k];
            Expr* v;
            if (wp_machine) {
                double d = mpfr_get_d(src, MPFR_RNDN);
                if (high && (k & 1)) d = -d;
                v = expr_new_real(d);
            } else {
                v = expr_new_mpfr_bits(wprec);
                mpfr_set(v->data.mpfr, src, MPFR_RNDN);
                if (high && (k & 1)) mpfr_neg(v->data.mpfr, v->data.mpfr, MPFR_RNDN);
            }
            rows[i] = wv_fn2(SYM_List, expr_new_integer(k), v);
        }
        r = wv_fn(SYM_List, rows, (size_t)taps);
        free(rows);
    }
    for (long k = 0; k < taps; k++) mpfr_clear(h[k]);
    free(h);
    return r;
}

/* DaubechiesWavelet[] = DaubechiesWavelet[2]. */
Expr* builtin_daubechieswavelet(Expr* res) {
    if (res->type != EXPR_FUNCTION || wv_len(res) != 0) return NULL;
    return wv_fn(SYM_DaubechiesWavelet, (Expr*[]){ expr_new_integer(2) }, 1);
}

#else  /* !USE_MPFR */
Expr* builtin_waveletphi(Expr* res) { (void)res; return NULL; }
Expr* builtin_waveletpsi(Expr* res) { (void)res; return NULL; }
Expr* builtin_waveletfiltercoefficients(Expr* res) { (void)res; return NULL; }
Expr* builtin_daubechieswavelet(Expr* res) { (void)res; return NULL; }
#endif

void wavelets_init(void) {
    symtab_add_builtin("WaveletPhi", builtin_waveletphi);
    symtab_add_builtin("WaveletPsi", builtin_waveletpsi);
    symtab_add_builtin("WaveletFilterCoefficients", builtin_waveletfiltercoefficients);
    symtab_add_builtin("DaubechiesWavelet", builtin_daubechieswavelet);
    static const char* heads[] = { "WaveletPhi", "WaveletPsi", "WaveletFilterCoefficients",
                                   "DaubechiesWavelet", "HaarWavelet" };
    for (size_t i = 0; i < sizeof heads / sizeof heads[0]; i++)
        symtab_get_def(heads[i])->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("WaveletPhi",
        "WaveletPhi[wave, x]\n"
        "\tthe scaling function phi(x) of the wavelet family wave (HaarWavelet[],\n"
        "\tDaubechiesWavelet[n]). Daubechies phi is the cascade's exact values on the\n"
        "\tdyadic grid of spacing 2^-MaxRecursion (default 8), linearly interpolated;\n"
        "\ta symbolic x gives Piecewise[{{InterpolatingFunction[..][x], 0 <= x <= 2n-1}}, 0].\n"
        "Options: MaxRecursion -> Automatic, WorkingPrecision -> MachinePrecision.");
    symtab_set_docstring("WaveletPsi",
        "WaveletPsi[wave, x]\n"
        "\tthe wavelet function psi(x) = 2 Sum_k g_k phi(2x - k), g_k = (-1)^k h_(1-k),\n"
        "\tof HaarWavelet[] or DaubechiesWavelet[n] (support [1-n, n]); as WaveletPhi.");
    symtab_set_docstring("WaveletFilterCoefficients",
        "WaveletFilterCoefficients[wave, \"PrimalLowpass\" | \"PrimalHighpass\"]\n"
        "\tthe filter {{k, h_k}, ...} of HaarWavelet[] or DaubechiesWavelet[n] (Sum h_k = 1),\n"
        "\tby spectral factorisation of the Daubechies polynomial. WorkingPrecision -> p.");
    symtab_set_docstring("DaubechiesWavelet",
        "DaubechiesWavelet[n] represents the Daubechies wavelet with n vanishing moments\n"
        "(2n filter taps); DaubechiesWavelet[] is DaubechiesWavelet[2].");
    symtab_set_docstring("HaarWavelet",
        "HaarWavelet[] represents the Haar wavelet: phi = 1 on [0, 1), psi = 1 on (0, 1/2)\n"
        "and -1 on (1/2, 1).");
}
