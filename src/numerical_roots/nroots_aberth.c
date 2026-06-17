/*
 * nroots_aberth.c — Aberth–Ehrlich simultaneous polynomial root-finder.
 *
 * The default NRoots engine.  All n roots are iterated together from Bini's
 * initial placement (starting points on circles whose radii come from the
 * upper convex hull / Newton polygon of (k, log|a_k|)), using the Aberth
 * correction
 *
 *     w_i = N_i / (1 - N_i * S_i),   N_i = p(z_i)/p'(z_i),
 *     S_i = sum_{j != i} 1/(z_i - z_j),   z_i <- z_i - w_i.
 *
 * Updates are applied in place (Aberth–Gauss-Seidel) for faster convergence.
 * Convergence is cubic at simple roots; roots of multiplicity m emerge as a
 * tight cluster (the caller snaps clusters to identical values).  References:
 * Ehrlich 1967; O. Aberth, Math. Comp. 27 (1973); D. Bini, Numer. Algorithms
 * 13 (1996) for the initialization.
 *
 * All arithmetic is MPFR complex (the ncpx toolkit).  Memory: every ncpx_init
 * pairs with an ncpx_clear; the caller owns the `roots` array.
 */

#include "nroots_internal.h"

#ifdef USE_MPFR

#include <math.h>
#include <stdlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ *
 *  Bini initial-point placement
 * ------------------------------------------------------------------ */
static void aberth_initial_points(const NrPoly* p, ncpx* roots, mpfr_prec_t wp) {
    int n = p->deg;

    /* log moduli l_k = log|c_k|, collected for nonzero coefficients only. */
    double* lk = (double*)malloc(sizeof(double) * (size_t)(n + 1));
    int*    kk = (int*)   malloc(sizeof(int)    * (size_t)(n + 1));
    int     np = 0;
    mpfr_t mag, lg;
    mpfr_init2(mag, wp); mpfr_init2(lg, wp);
    for (int k = 0; k <= n; k++) {
        ncpx_abs(mag, &p->c[k]);
        if (mpfr_zero_p(mag)) continue;        /* skip zero coefficients */
        mpfr_log(lg, mag, MPFR_RNDN);
        lk[np] = mpfr_get_d(lg, MPFR_RNDN);
        kk[np] = k;
        np++;
    }
    mpfr_clear(mag); mpfr_clear(lg);

    /* Upper convex hull of (kk[i], lk[i]) via monotone chain (pop while the
     * turn is non-clockwise: cross >= 0). */
    int* hull = (int*)malloc(sizeof(int) * (size_t)np);
    int  hn = 0;
    for (int i = 0; i < np; i++) {
        while (hn >= 2) {
            int a = hull[hn - 2], b = hull[hn - 1], c = i;
            double cross = (double)(kk[b] - kk[a]) * (lk[c] - lk[a])
                         - (lk[b] - lk[a]) * (double)(kk[c] - kk[a]);
            if (cross >= 0.0) hn--; else break;
        }
        hull[hn++] = i;
    }

    /* Walk hull edges; each edge of horizontal length dk seeds dk points on a
     * circle of radius exp((l_a - l_b)/dk). */
    int idx = 0;
    double sigma = 0.7;   /* global rotation, keeps starts off the real axis */
    for (int e = 0; e + 1 < hn; e++) {
        int ia = hull[e], ib = hull[e + 1];
        int dk = kk[ib] - kk[ia];
        if (dk <= 0) continue;
        double r = exp((lk[ia] - lk[ib]) / (double)dk);
        double edgephase = 2.0 * M_PI * (double)e / (double)(n > 0 ? n : 1);
        for (int t = 0; t < dk && idx < n; t++) {
            double theta = (2.0 * M_PI / (double)dk) * (double)t + sigma + edgephase;
            ncpx_set_d(&roots[idx], r * cos(theta), r * sin(theta));
            idx++;
        }
    }
    /* Defensive: if rounding left points unfilled, scatter them on a unit
     * circle (should not happen for a deflated polynomial). */
    for (; idx < n; idx++) {
        double theta = 2.0 * M_PI * (double)idx / (double)n + sigma;
        ncpx_set_d(&roots[idx], cos(theta), sin(theta));
    }

    free(lk); free(kk); free(hull);
}

/* ------------------------------------------------------------------ *
 *  Aberth iteration
 * ------------------------------------------------------------------ */
int nr_aberth(const NrPoly* p, ncpx* roots, int max_iter, mpfr_prec_t wp) {
    int n = p->deg;

    /* Linear case. */
    if (n == 1) {
        ncpx_div(&roots[0], &p->c[0], &p->c[1], wp);  /* -c0/c1 */
        mpfr_neg(roots[0].re, roots[0].re, MPFR_RNDN);
        mpfr_neg(roots[0].im, roots[0].im, MPFR_RNDN);
        return 0;
    }

    aberth_initial_points(p, roots, wp);

    int* done = (int*)calloc((size_t)n, sizeof(int));
    ncpx newt, S, diff, inv, prod, denom, w, one;
    ncpx val, der;
    ncpx_init(&newt, wp); ncpx_init(&S, wp); ncpx_init(&diff, wp);
    ncpx_init(&inv, wp); ncpx_init(&prod, wp); ncpx_init(&denom, wp);
    ncpx_init(&w, wp); ncpx_init(&one, wp);
    ncpx_init(&val, wp); ncpx_init(&der, wp);
    ncpx_set_ui(&one, 1);

    mpfr_t tol, wabs, zabs, bound, dmag;
    mpfr_init2(tol, wp); mpfr_init2(wabs, wp); mpfr_init2(zabs, wp);
    mpfr_init2(bound, wp); mpfr_init2(dmag, wp);
    mpfr_set_ui(tol, 1, MPFR_RNDN);
    mpfr_div_2si(tol, tol, (long)wp - 8, MPFR_RNDN);   /* 2^-(wp-8) */

    int converged = 0;
    for (int it = 0; it < max_iter && !converged; it++) {
        int all_done = 1;
        for (int i = 0; i < n; i++) {
            if (done[i]) continue;
            all_done = 0;

            nr_poly_eval(p, &roots[i], &val, &der, wp);
            ncpx_abs(dmag, &der);
            if (mpfr_zero_p(dmag)) {
                /* Stationary point: nudge and retry next sweep. */
                mpfr_add_d(roots[i].re, roots[i].re, 1e-3, MPFR_RNDN);
                continue;
            }
            ncpx_div(&newt, &val, &der, wp);          /* N_i = p/p' */

            /* S_i = sum_{j != i} 1/(z_i - z_j). */
            ncpx_set_ui(&S, 0);
            for (int j = 0; j < n; j++) {
                if (j == i) continue;
                ncpx_sub(&diff, &roots[i], &roots[j]);
                ncpx_abs(bound, &diff);
                if (mpfr_zero_p(bound)) continue;     /* coincident; skip term */
                ncpx_div(&inv, &one, &diff, wp);
                ncpx_add(&S, &S, &inv);
            }

            /* w_i = N_i / (1 - N_i * S_i). */
            ncpx_mul(&prod, &newt, &S, wp);
            ncpx_sub(&denom, &one, &prod);
            ncpx_abs(bound, &denom);
            if (mpfr_zero_p(bound)) ncpx_set(&w, &newt);
            else ncpx_div(&w, &newt, &denom, wp);

            ncpx_sub(&roots[i], &roots[i], &w);

            ncpx_abs(wabs, &w);
            ncpx_abs(zabs, &roots[i]);
            mpfr_set_ui(bound, 1, MPFR_RNDN);
            if (mpfr_cmp(zabs, bound) > 0) mpfr_set(bound, zabs, MPFR_RNDN);
            mpfr_mul(bound, bound, tol, MPFR_RNDN);
            if (mpfr_cmp(wabs, bound) <= 0) done[i] = 1;
        }
        if (all_done) converged = 1;
    }

    /* Final polish on the original polynomial (cleans up any laggards). */
    for (int i = 0; i < n; i++) nr_newton_polish(p, &roots[i], wp, 8);

    ncpx_clear(&newt); ncpx_clear(&S); ncpx_clear(&diff);
    ncpx_clear(&inv); ncpx_clear(&prod); ncpx_clear(&denom);
    ncpx_clear(&w); ncpx_clear(&one);
    ncpx_clear(&val); ncpx_clear(&der);
    mpfr_clear(tol); mpfr_clear(wabs); mpfr_clear(zabs);
    mpfr_clear(bound); mpfr_clear(dmag);
    free(done);
    return 0;
}

#endif /* USE_MPFR */
