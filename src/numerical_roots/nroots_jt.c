/*
 * nroots_jt.c — Jenkins–Traub polynomial root-finder.
 *
 * Implements the three-stage shifted-deflation algorithm of Jenkins & Traub.
 * The complex variant (CPOLY, ACM TOMS Algorithm 419, CACM 1972) is used for
 * all coefficient types: it natively handles complex coefficients, and on a
 * real-coefficient polynomial it returns the same roots as the real-arithmetic
 * specialization (RPOLY, ACM TOMS 493) — the latter is only a speed/storage
 * optimization, so it is folded into the complex path here.
 *
 * One root is extracted at a time:
 *   Stage 1 (no-shift):    H <- (H - (H(0)/P(0)) P) / z, a few times, to
 *                          emphasize the smallest-modulus zero.
 *   Stage 2 (fixed-shift): H is shaped with a fixed shift s near a zero
 *                          (modulus = Cauchy lower bound, sweeping angles).
 *   Stage 3 (variable):    s_{k+1} = s_k - P(s_k)/H(s_k), cubically convergent.
 * The found zero is polished against the ORIGINAL polynomial, recorded, then
 * deflated out; a zero of multiplicity m is found m times.
 *
 * All arithmetic is MPFR complex (ncpx).  Memory: every ncpx_init pairs with an
 * ncpx_clear; the caller owns the `roots` array.
 */

#include "nroots_internal.h"

#ifdef USE_MPFR

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- small ncpx-array helpers ------------------------------------- */

/* out = c(s), Horner over a degree-`deg` complex coefficient array. */
static void jt_eval(const ncpx* c, int deg, const ncpx* s, ncpx* out, mpfr_prec_t wp) {
    ncpx b; ncpx_init(&b, wp);
    ncpx_set(&b, &c[deg]);
    for (int i = deg - 1; i >= 0; i--) { ncpx_mul(&b, &b, s, wp); ncpx_add(&b, &b, &c[i]); }
    ncpx_set(out, &b);
    ncpx_clear(&b);
}

/* H <- (H - (H(s)/P(s)) P) / (z - s), the Jenkins–Traub H recurrence with
 * shift s.  P has degree n, H degree n-1.  Returns 0 if P(s) ~ 0 (skip). */
static int jt_next_h(const ncpx* P, int n, ncpx* H, const ncpx* s, mpfr_prec_t wp) {
    ncpx Ps, Hs, c;
    ncpx_init(&Ps, wp); ncpx_init(&Hs, wp); ncpx_init(&c, wp);
    jt_eval(P, n, s, &Ps, wp);
    jt_eval(H, n - 1, s, &Hs, wp);
    mpfr_t m; mpfr_init2(m, wp); ncpx_abs(m, &Ps);
    if (mpfr_zero_p(m)) { mpfr_clear(m); ncpx_clear(&Ps); ncpx_clear(&Hs); ncpx_clear(&c); return 0; }
    ncpx_div(&c, &Hs, &Ps, wp);                 /* c = H(s)/P(s) */

    /* N = H - c P  (degree n: N[n] = -c P[n], N[i] = H[i]-c P[i] for i<n). */
    ncpx* N = (ncpx*)malloc(sizeof(ncpx) * (size_t)(n + 1));
    ncpx t; ncpx_init(&t, wp);
    for (int i = 0; i < n; i++) {
        ncpx_init(&N[i], wp);
        ncpx_mul(&t, &c, &P[i], wp);
        ncpx_sub(&N[i], &H[i], &t);
    }
    ncpx_init(&N[n], wp);
    ncpx_mul(&N[n], &c, &P[n], wp);
    ncpx_neg(&N[n], &N[n]);

    /* H <- N / (z - s) via synthetic division (quotient degree n-1). */
    ncpx q; ncpx_init(&q, wp);
    ncpx_set(&q, &N[n]);                          /* q[n-1] */
    ncpx_set(&H[n - 1], &q);
    for (int i = n - 1; i >= 1; i--) {
        ncpx_mul(&t, &q, s, wp);
        ncpx_add(&q, &N[i], &t);                  /* q[i-1] = N[i] + s*q[i] */
        ncpx_set(&H[i - 1], &q);
    }

    /* Normalize H to monic.  The recurrence is linear in H, so scaling leaves
     * its direction unchanged; keeping it monic bounds magnitudes and makes the
     * variable-shift iterate  s - P(s)/H(s)  correctly scaled. */
    ncpx lead; ncpx_init(&lead, wp);
    ncpx_set(&lead, &H[n - 1]);
    ncpx_abs(m, &lead);
    if (!mpfr_zero_p(m)) {
        for (int i = 0; i < n; i++) ncpx_div(&H[i], &H[i], &lead, wp);
    }
    ncpx_clear(&lead);

    for (int i = 0; i <= n; i++) ncpx_clear(&N[i]);
    free(N);
    ncpx_clear(&t); ncpx_clear(&q);
    mpfr_clear(m); ncpx_clear(&Ps); ncpx_clear(&Hs); ncpx_clear(&c);
    return 1;
}

/* Cauchy lower bound on the smallest root modulus: the unique positive root of
 * x^n + |a_{n-1}| x^{n-1} + ... + |a_1| x - |a_0| = 0.  Computed in double
 * precision (it only sets the shift radius). */
static double jt_cauchy_radius(const ncpx* P, int n, mpfr_prec_t wp) {
    double* a = (double*)malloc(sizeof(double) * (size_t)(n + 1));
    mpfr_t m; mpfr_init2(m, wp);
    for (int i = 0; i <= n; i++) { ncpx_abs(m, &P[i]); a[i] = mpfr_get_d(m, MPFR_RNDN); }
    mpfr_clear(m);
    /* g(x) = -a0 + a1 x + ... + a_{n-1} x^{n-1} + a_n x^n. */
    double an = a[n] > 0 ? a[n] : 1.0;
    double x = (a[0] > 0) ? pow(a[0] / an, 1.0 / (double)n) : 1.0;
    if (!(x > 0.0)) x = 1.0;
    for (int it = 0; it < 60; it++) {
        double g = a[n], dg = 0.0;            /* Horner with negated a0 */
        for (int i = n - 1; i >= 0; i--) {
            double ci = (i == 0) ? -a[0] : a[i];
            dg = dg * x + g;
            g = g * x + ci;
        }
        if (dg == 0.0) break;
        double xn = x - g / dg;
        if (xn <= 0.0) xn = x * 0.5;
        if (fabs(xn - x) < 1e-4 * x) { x = xn; break; }
        x = xn;
    }
    free(a);
    if (!(x > 0.0) || !isfinite(x)) x = 1.0;
    return x;
}

/* Attempt to find one zero of P (degree n) starting from fixed shift s0.
 * On success writes the zero to *zero and returns 1. */
static int jt_one_zero(const ncpx* P, int n, ncpx* H, const ncpx* s0,
                       ncpx* zero, mpfr_prec_t wp) {
    ncpx s, Ps, Hs, ds;
    ncpx_init(&s, wp); ncpx_init(&Ps, wp); ncpx_init(&Hs, wp); ncpx_init(&ds, wp);
    mpfr_t step, sc, tol, hm; mpfr_init2(step, wp); mpfr_init2(sc, wp);
    mpfr_init2(tol, wp); mpfr_init2(hm, wp);
    mpfr_set_ui(tol, 1, MPFR_RNDN);
    mpfr_div_2si(tol, tol, (long)wp - 6, MPFR_RNDN);

    /* Residual gate (relative to the coefficient scale). */
    mpfr_t maxc, resid, rtol; mpfr_init2(maxc, wp); mpfr_init2(resid, wp); mpfr_init2(rtol, wp);
    mpfr_set_zero(maxc, 1);
    for (int i = 0; i <= n; i++) { ncpx_abs(hm, &P[i]); if (mpfr_cmp(hm, maxc) > 0) mpfr_set(maxc, hm, MPFR_RNDN); }
    mpfr_set(rtol, maxc, MPFR_RNDN);
    mpfr_div_2si(rtol, rtol, (long)wp - 12, MPFR_RNDN);

    /* Stage 2: fixed shift s0 (shape H). */
    for (int i = 0; i < 10; i++) jt_next_h(P, n, H, s0, wp);

    /* Stage 3: variable shift, cubically convergent. */
    ncpx_set(&s, s0);
    int ok = 0;
    for (int it = 0; it < 200; it++) {
        jt_next_h(P, n, H, &s, wp);
        jt_eval(P, n, &s, &Ps, wp);
        jt_eval(H, n - 1, &s, &Hs, wp);
        ncpx_abs(hm, &Hs);
        if (mpfr_zero_p(hm) || !mpfr_number_p(Hs.re) || !mpfr_number_p(Hs.im)) break;
        ncpx_div(&ds, &Ps, &Hs, wp);
        if (!mpfr_number_p(ds.re) || !mpfr_number_p(ds.im)) break;
        ncpx_sub(&s, &s, &ds);
        if (!mpfr_number_p(s.re) || !mpfr_number_p(s.im)) break;

        ncpx_abs(step, &ds);
        ncpx_abs(sc, &s);
        mpfr_t one; mpfr_init2(one, wp); mpfr_set_ui(one, 1, MPFR_RNDN);
        if (mpfr_cmp(sc, one) < 0) mpfr_set(sc, one, MPFR_RNDN);
        mpfr_mul(sc, sc, tol, MPFR_RNDN);
        int small = mpfr_cmp(step, sc) <= 0;
        mpfr_clear(one);
        if (small) {
            jt_eval(P, n, &s, &Ps, wp);
            ncpx_abs(resid, &Ps);
            if (mpfr_cmp(resid, rtol) <= 0) ok = 1;
            break;
        }
    }
    if (ok) ncpx_set(zero, &s);
    mpfr_clear(maxc); mpfr_clear(resid); mpfr_clear(rtol);

    ncpx_clear(&s); ncpx_clear(&Ps); ncpx_clear(&Hs); ncpx_clear(&ds);
    mpfr_clear(step); mpfr_clear(sc); mpfr_clear(tol); mpfr_clear(hm);
    return ok;
}

/* Deflate P (degree n, monic not required) by (z - r): quotient degree n-1
 * via synthetic division; writes into Pout (length n). */
static void jt_deflate(const ncpx* P, int n, const ncpx* r, ncpx* Pout, mpfr_prec_t wp) {
    ncpx q, t; ncpx_init(&q, wp); ncpx_init(&t, wp);
    ncpx_set(&q, &P[n]);
    ncpx_set(&Pout[n - 1], &q);
    for (int i = n - 1; i >= 1; i--) {
        ncpx_mul(&t, &q, r, wp);
        ncpx_add(&q, &P[i], &t);
        ncpx_set(&Pout[i - 1], &q);
    }
    ncpx_clear(&q); ncpx_clear(&t);
}

int nr_jenkinstraub(const NrPoly* p, ncpx* roots, int max_iter, mpfr_prec_t wp) {
    (void)max_iter;
    int n0 = p->deg;

    /* Working (deflating) copy of the coefficients, made monic — the Jenkins–
     * Traub iterate s - P(s)/H(s) assumes a monic P (roots are unchanged). */
    ncpx* P = (ncpx*)malloc(sizeof(ncpx) * (size_t)(n0 + 1));
    for (int i = 0; i <= n0; i++) { ncpx_init(&P[i], wp); ncpx_set(&P[i], &p->c[i]); }
    {
        ncpx lead; ncpx_init(&lead, wp); ncpx_set(&lead, &P[n0]);
        for (int i = 0; i <= n0; i++) ncpx_div(&P[i], &P[i], &lead, wp);
        ncpx_clear(&lead);
    }

    ncpx* H = (ncpx*)malloc(sizeof(ncpx) * (size_t)(n0));   /* degree n-1 max */
    for (int i = 0; i < n0; i++) ncpx_init(&H[i], wp);

    ncpx s0, zero, czero;
    ncpx_init(&s0, wp); ncpx_init(&zero, wp); ncpx_init(&czero, wp);
    ncpx_set_d(&zero, 0.0, 0.0);     /* ncpx_init leaves NaN; make these 0 */
    ncpx_set_d(&czero, 0.0, 0.0);    /* the fixed no-shift origin */

    int found = 0;
    int n = n0;
    int rc = 0;
    while (n >= 1 && rc == 0) {
        if (n == 1) {
            /* root = -P0/P1 */
            ncpx_div(&roots[found], &P[0], &P[1], wp);
            mpfr_neg(roots[found].re, roots[found].re, MPFR_RNDN);
            mpfr_neg(roots[found].im, roots[found].im, MPFR_RNDN);
            nr_newton_polish(p, &roots[found], wp, 12);
            found++;
            break;
        }

        /* Stage 1: no-shift, build H from P'/n. */
        for (int i = 0; i < n; i++) {
            /* H[i] = (i+1) P[i+1] / n */
            ncpx_set(&H[i], &P[i + 1]);
            mpfr_mul_ui(H[i].re, H[i].re, (unsigned long)(i + 1), MPFR_RNDN);
            mpfr_mul_ui(H[i].im, H[i].im, (unsigned long)(i + 1), MPFR_RNDN);
            mpfr_div_ui(H[i].re, H[i].re, (unsigned long)n, MPFR_RNDN);
            mpfr_div_ui(H[i].im, H[i].im, (unsigned long)n, MPFR_RNDN);
        }
        for (int i = 0; i < 5; i++) jt_next_h(P, n, H, &czero, wp);

        /* Stage 2/3: sweep shift angles on the Cauchy circle. */
        double beta = jt_cauchy_radius(P, n, wp);
        int got = 0;
        double angle = 49.0 * M_PI / 180.0;
        for (int tries = 0; tries < 12 && !got; tries++) {
            ncpx_set_d(&s0, beta * cos(angle), beta * sin(angle));
            /* Rebuild H from P' for each fresh shift attempt. */
            for (int i = 0; i < n; i++) {
                ncpx_set(&H[i], &P[i + 1]);
                mpfr_mul_ui(H[i].re, H[i].re, (unsigned long)(i + 1), MPFR_RNDN);
                mpfr_mul_ui(H[i].im, H[i].im, (unsigned long)(i + 1), MPFR_RNDN);
                mpfr_div_ui(H[i].re, H[i].re, (unsigned long)n, MPFR_RNDN);
                mpfr_div_ui(H[i].im, H[i].im, (unsigned long)n, MPFR_RNDN);
            }
            for (int i = 0; i < 5; i++) jt_next_h(P, n, H, &czero, wp);
            if (jt_one_zero(P, n, H, &s0, &zero, wp)) got = 1;
            angle += 94.0 * M_PI / 180.0;
        }
        if (!got) { rc = -1; break; }

        /* Polish against the original polynomial, record, deflate. */
        nr_newton_polish(p, &zero, wp, 16);
        ncpx_set(&roots[found], &zero);
        found++;

        ncpx* Pn = (ncpx*)malloc(sizeof(ncpx) * (size_t)n);   /* degree n-1 */
        for (int i = 0; i < n; i++) ncpx_init(&Pn[i], wp);
        jt_deflate(P, n, &zero, Pn, wp);
        for (int i = 0; i <= n; i++) ncpx_clear(&P[i]);
        free(P);
        P = Pn;
        n--;
    }

    if (rc == 0 && found != n0) rc = -1;

    for (int i = 0; i <= n; i++) ncpx_clear(&P[i]);   /* P has degree n now */
    free(P);
    for (int i = 0; i < n0; i++) ncpx_clear(&H[i]);
    free(H);
    ncpx_clear(&s0); ncpx_clear(&zero); ncpx_clear(&czero);
    return rc;
}

#endif /* USE_MPFR */
