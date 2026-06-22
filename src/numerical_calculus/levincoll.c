/*
 * levincoll.c — Levin collocation rule for highly oscillatory integrals.
 *
 * See levincoll.h for the mathematics.  This file implements the machine-
 * precision 1-D core: build the Chebyshev collocation system for the Levin ODE
 * p' + i g' p = f at Chebyshev–Gauss–Lobatto nodes, solve it with LAPACK's
 * complex LU, and evaluate the boundary term.  The order is refined (doubled)
 * until two successive estimates agree.
 */
#include "levincoll.h"

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "../linalg/lapack.h"
#ifdef USE_MPFR
#include "../numeric_complex.h"
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Smallest collocation order tried, then doubled. */
#define LEVIN_N_MIN 4
/* Reciprocal-condition floor: below this the collocation matrix is effectively
 * singular (weak oscillation / stationary phase) and Levin does not apply. */
#define LEVIN_RCOND_MIN 1e-13

/* Read the (row r, col c) entry of a column-major interleaved-complex matrix. */
static inline double _Complex zget(const double* A, int n, int r, int c) {
    const double* p = A + 2 * ((size_t)c * n + r);
    return p[0] + p[1] * I;
}

/* Chebyshev T_k(t) and dT_k/dt for k = 0..n-1 at a single t in [-1,1].
 * T_0=1, T_1=t, T_k = 2t T_{k-1} - T_{k-2};  T_k'(t) = k U_{k-1}(t),
 * U_0=1, U_1=2t, U_k = 2t U_{k-1} - U_{k-2}.  At the endpoints t=±1 the
 * derivative recurrence is numerically delicate, so use the closed forms
 * T_k'(1) = k^2 and T_k'(-1) = (-1)^{k+1} k^2 there. */
static void levin_cheb_row(double t, int n, double* T, double* dT) {
    T[0] = 1.0;
    if (n > 1) T[1] = t;
    for (int k = 2; k < n; k++) T[k] = 2.0 * t * T[k - 1] - T[k - 2];

    bool at_hi = (t >= 1.0), at_lo = (t <= -1.0);
    if (at_hi || at_lo) {
        for (int k = 0; k < n; k++) {
            double k2 = (double)k * (double)k;
            dT[k] = at_hi ? k2 : ((k % 2 == 0) ? -k2 : k2);
        }
        return;
    }
    /* Interior: U_{k-1}(t), then T_k'(t) = k U_{k-1}(t). */
    double Ukm1 = 0.0;       /* U_{-1} = 0 */
    double Uk = 1.0;         /* U_0 = 1   */
    dT[0] = 0.0;
    for (int k = 1; k < n; k++) {
        dT[k] = (double)k * Uk;   /* U_{k-1} */
        double Ukp1 = 2.0 * t * Uk - Ukm1;
        Ukm1 = Uk; Uk = Ukp1;
    }
}

/* Solve A x = b in place (b -> x) given the LU factors produced by zgetrf in
 * the column-major interleaved buffer A and the 1-indexed pivots ipiv.  A holds
 * unit-lower L (strict lower triangle) and upper U (incl. diagonal).  Returns
 * false if a U diagonal underflows (singular). */
static bool levin_zsolve_lu(int n, const double* A, const int* ipiv,
                            double _Complex* b) {
    /* Apply the row interchanges P from the factorisation: for i ascending,
     * swap b[i] with b[ipiv[i]-1]. */
    for (int i = 0; i < n; i++) {
        int p = ipiv[i] - 1;
        if (p != i) { double _Complex t = b[i]; b[i] = b[p]; b[p] = t; }
    }
    /* Forward substitution: L y = Pb, L unit lower triangular. */
    for (int i = 1; i < n; i++) {
        double _Complex s = b[i];
        for (int k = 0; k < i; k++) s -= zget(A, n, i, k) * b[k];
        b[i] = s;
    }
    /* Back substitution: U x = y. */
    for (int i = n - 1; i >= 0; i--) {
        double _Complex s = b[i];
        for (int k = i + 1; k < n; k++) s -= zget(A, n, i, k) * b[k];
        double _Complex d = zget(A, n, i, i);
        if (cabs(d) < 1e-300) return false;
        b[i] = s / d;
    }
    return true;
}

/* One Levin solve at collocation order n.  Writes the kernel-projected integral
 * estimate to *out and the matrix's reciprocal condition number to *rcond.
 * Returns false on a sampling failure, a singular system, or no LAPACK. */
static bool levin_solve_order(int n, double a, double b,
                              GkSampleMachine amp,    void* amp_ctx,
                              GkSampleMachine gprime, void* gprime_ctx,
                              GkSampleMachine gphase, void* gphase_ctx,
                              LevinKernel kernel,
                              double _Complex* out, double* rcond) {
    const double half = 0.5 * (b - a);
    const double mid  = 0.5 * (a + b);
    const double dtdx = 1.0 / half;        /* dt/dx for t = (x-mid)/half */

    double* A   = malloc(sizeof(double) * 2 * (size_t)n * n);
    double _Complex* rhs = malloc(sizeof(double _Complex) * n);
    double* T   = malloc(sizeof(double) * n);
    double* dT  = malloc(sizeof(double) * n);
    int*    ipiv = malloc(sizeof(int) * n);
    bool ok = (A && rhs && T && dT && ipiv);

    for (int j = 0; ok && j < n; j++) {
        double t = cos(M_PI * (double)j / (double)(n - 1));   /* CGL node */
        double x = mid + half * t;
        double _Complex fj, gpj;
        if (!amp(amp_ctx, x, &fj) || !gprime(gprime_ctx, x, &gpj)) { ok = false; break; }
        levin_cheb_row(t, n, T, dT);
        for (int k = 0; k < n; k++) {
            /* A[j,k] = u_k'(x) + i g'(x) u_k(x),  u_k(x) = T_k(t(x)). */
            double _Complex v = (dT[k] * dtdx) + I * gpj * T[k];
            double* p = A + 2 * ((size_t)k * n + j);
            p[0] = creal(v); p[1] = cimag(v);
        }
        rhs[j] = fj;
    }

    bool solved = false;
    double _Complex Ival = 0.0;
    if (ok) {
        double anorm = mat_lapack_zlange('1', n, n, A, n);
        int info = mat_lapack_zgetrf(n, n, A, n, ipiv);
        if (info == 0 && anorm >= 0.0) {
            double rc = 0.0;
            int ci = mat_lapack_zgecon('1', n, A, n, anorm, &rc);
            *rcond = (ci == 0) ? rc : 0.0;
            if (*rcond >= LEVIN_RCOND_MIN && levin_zsolve_lu(n, A, ipiv, rhs)) {
                /* p(b): t=+1 so T_k(1)=1; p(a): t=-1 so T_k(-1)=(-1)^k. */
                double _Complex pb = 0.0, pa = 0.0;
                for (int k = 0; k < n; k++) {
                    pb += rhs[k];
                    pa += (k % 2 == 0) ? rhs[k] : -rhs[k];
                }
                double _Complex ga, gb;
                if (gphase(gphase_ctx, a, &ga) && gphase(gphase_ctx, b, &gb)) {
                    Ival = pb * cexp(I * gb) - pa * cexp(I * ga);
                    solved = true;
                }
            }
        } else {
            *rcond = 0.0;
        }
    }

    free(A); free(rhs); free(T); free(dT); free(ipiv);
    if (!solved) return false;

    switch (kernel) {
        case LEVIN_KERNEL_COS: *out = creal(Ival); break;
        case LEVIN_KERNEL_SIN: *out = cimag(Ival); break;
        default:               *out = Ival;        break;
    }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Prepared (factored-once) solver                                    *
 * ------------------------------------------------------------------ */

struct LevinPrep {
    int     n;
    double  a, b;
    double* A;       /* order-n LU (col-major interleaved), 2*n*n doubles */
    int*    ipiv;
    double* nodes;   /* the n collocation abscissae x_j                   */
};

LevinPrep* levin_prepare_machine(double a, double b, int n,
                                 GkSampleMachine gprime, void* gprime_ctx) {
    if (n < 2 || !(b > a)) return NULL;
    const double half = 0.5 * (b - a);
    const double mid  = 0.5 * (a + b);
    const double dtdx = 1.0 / half;

    LevinPrep* P = malloc(sizeof(LevinPrep));
    double* A      = malloc(sizeof(double) * 2 * (size_t)n * n);
    int*    ipiv   = malloc(sizeof(int) * n);
    double* nodes  = malloc(sizeof(double) * n);
    double* T      = malloc(sizeof(double) * n);
    double* dT     = malloc(sizeof(double) * n);
    bool ok = (P && A && ipiv && nodes && T && dT);

    for (int j = 0; ok && j < n; j++) {
        double t = cos(M_PI * (double)j / (double)(n - 1));
        double x = mid + half * t;
        nodes[j] = x;
        double _Complex gpj;
        if (!gprime(gprime_ctx, x, &gpj)) { ok = false; break; }
        levin_cheb_row(t, n, T, dT);
        for (int k = 0; k < n; k++) {
            double _Complex v = (dT[k] * dtdx) + I * gpj * T[k];
            double* p = A + 2 * ((size_t)k * n + j);
            p[0] = creal(v); p[1] = cimag(v);
        }
    }
    if (ok) {
        double anorm = mat_lapack_zlange('1', n, n, A, n);
        int info = mat_lapack_zgetrf(n, n, A, n, ipiv);
        if (info != 0 || anorm < 0.0) ok = false;
        if (ok) {
            double rc = 0.0;
            int ci = mat_lapack_zgecon('1', n, A, n, anorm, &rc);
            if (ci != 0 || rc < LEVIN_RCOND_MIN) ok = false;
        }
    }
    free(T); free(dT);
    if (!ok) { free(A); free(ipiv); free(nodes); free(P); return NULL; }
    P->n = n; P->a = a; P->b = b; P->A = A; P->ipiv = ipiv; P->nodes = nodes;
    return P;
}

bool levin_prepared_solve(const LevinPrep* P,
                          GkSampleMachine amp,    void* amp_ctx,
                          GkSampleMachine gphase, void* gphase_ctx,
                          LevinKernel kernel, double _Complex* out) {
    int n = P->n;
    double _Complex* rhs = malloc(sizeof(double _Complex) * n);
    if (!rhs) return false;
    for (int j = 0; j < n; j++)
        if (!amp(amp_ctx, P->nodes[j], &rhs[j])) { free(rhs); return false; }
    if (!levin_zsolve_lu(n, P->A, P->ipiv, rhs)) { free(rhs); return false; }

    double _Complex pb = 0.0, pa = 0.0;
    for (int k = 0; k < n; k++) {
        pb += rhs[k];
        pa += (k % 2 == 0) ? rhs[k] : -rhs[k];
    }
    free(rhs);
    double _Complex ga, gb;
    if (!gphase(gphase_ctx, P->a, &ga) || !gphase(gphase_ctx, P->b, &gb)) return false;
    double _Complex Ival = pb * cexp(I * gb) - pa * cexp(I * ga);
    if (!isfinite(creal(Ival)) || !isfinite(cimag(Ival))) return false;
    switch (kernel) {
        case LEVIN_KERNEL_COS: *out = creal(Ival); break;
        case LEVIN_KERNEL_SIN: *out = cimag(Ival); break;
        default:               *out = Ival;        break;
    }
    return true;
}

void levin_prepare_free(LevinPrep* P) {
    if (!P) return;
    free(P->A); free(P->ipiv); free(P->nodes); free(P);
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------ *
 *  Arbitrary-precision path (ncpx complex, in-house LU)               *
 * ------------------------------------------------------------------ */

/* Chebyshev T_k(t) and dT_k/dt at one t in [-1,1], all at precision wp.  Uses
 * the same recurrences as the machine path, with the endpoint closed forms
 * T_k'(±1) for numerical cleanliness. */
static void levin_cheb_row_mpfr(const mpfr_t t, int n, mpfr_t* T, mpfr_t* dT,
                                mpfr_prec_t wp) {
    mpfr_set_ui(T[0], 1, MPFR_RNDN);
    if (n > 1) mpfr_set(T[1], t, MPFR_RNDN);
    mpfr_t tmp; mpfr_init2(tmp, wp);
    for (int k = 2; k < n; k++) {            /* T_k = 2 t T_{k-1} - T_{k-2} */
        mpfr_mul(tmp, t, T[k - 1], MPFR_RNDN);
        mpfr_mul_2ui(tmp, tmp, 1, MPFR_RNDN);
        mpfr_sub(T[k], tmp, T[k - 2], MPFR_RNDN);
    }
    int at_hi = (mpfr_cmp_si(t, 1) >= 0), at_lo = (mpfr_cmp_si(t, -1) <= 0);
    if (at_hi || at_lo) {
        for (int k = 0; k < n; k++) {
            mpfr_set_si(dT[k], (long)k * (long)k, MPFR_RNDN);
            if (at_lo && (k % 2 == 0)) mpfr_neg(dT[k], dT[k], MPFR_RNDN);
        }
        mpfr_clear(tmp); return;
    }
    /* Interior: U_{k-1}(t), then T_k'(t) = k U_{k-1}(t). */
    mpfr_t Ukm1, Uk, Ukp1; mpfr_init2(Ukm1, wp); mpfr_init2(Uk, wp); mpfr_init2(Ukp1, wp);
    mpfr_set_ui(Ukm1, 0, MPFR_RNDN);    /* U_{-1} */
    mpfr_set_ui(Uk, 1, MPFR_RNDN);      /* U_0    */
    mpfr_set_ui(dT[0], 0, MPFR_RNDN);
    for (int k = 1; k < n; k++) {
        mpfr_mul_ui(dT[k], Uk, (unsigned long)k, MPFR_RNDN);
        mpfr_mul(Ukp1, t, Uk, MPFR_RNDN);
        mpfr_mul_2ui(Ukp1, Ukp1, 1, MPFR_RNDN);
        mpfr_sub(Ukp1, Ukp1, Ukm1, MPFR_RNDN);
        mpfr_set(Ukm1, Uk, MPFR_RNDN); mpfr_set(Uk, Ukp1, MPFR_RNDN);
    }
    mpfr_clear(tmp); mpfr_clear(Ukm1); mpfr_clear(Uk); mpfr_clear(Ukp1);
}

/* Solve A c = b for an n×n ncpx system (row-major) by Gaussian elimination
 * with partial pivoting, in place (b becomes c).  Returns false if singular. */
static bool levin_ncpx_solve(ncpx* A, ncpx* b, int n, mpfr_prec_t wp) {
    mpfr_t mag, best; mpfr_init2(mag, wp); mpfr_init2(best, wp);
    ncpx factor, tmp; ncpx_init(&factor, wp); ncpx_init(&tmp, wp);
    bool ok = true;
    for (int col = 0; col < n && ok; col++) {
        int piv = col; mpfr_set_ui(best, 0, MPFR_RNDN);
        for (int r = col; r < n; r++) {
            ncpx_abs(mag, &A[r * n + col]);
            if (mpfr_cmp(mag, best) > 0) { mpfr_set(best, mag, MPFR_RNDN); piv = r; }
        }
        if (mpfr_zero_p(best)) { ok = false; break; }
        if (piv != col) {
            for (int k = 0; k < n; k++) {
                mpfr_swap(A[piv * n + k].re, A[col * n + k].re);
                mpfr_swap(A[piv * n + k].im, A[col * n + k].im);
            }
            mpfr_swap(b[piv].re, b[col].re); mpfr_swap(b[piv].im, b[col].im);
        }
        for (int r = col + 1; r < n; r++) {
            ncpx_div(&factor, &A[r * n + col], &A[col * n + col], wp);
            for (int k = col; k < n; k++) {
                ncpx_mul(&tmp, &factor, &A[col * n + k], wp);
                ncpx_sub(&A[r * n + k], &A[r * n + k], &tmp);
            }
            ncpx_mul(&tmp, &factor, &b[col], wp);
            ncpx_sub(&b[r], &b[r], &tmp);
        }
    }
    for (int i = n - 1; i >= 0 && ok; i--) {
        for (int k = i + 1; k < n; k++) {
            ncpx_mul(&tmp, &A[i * n + k], &b[k], wp);
            ncpx_sub(&b[i], &b[i], &tmp);
        }
        ncpx_div(&b[i], &b[i], &A[i * n + i], wp);
    }
    mpfr_clear(mag); mpfr_clear(best); ncpx_clear(&factor); ncpx_clear(&tmp);
    return ok;
}

/* One arbitrary-precision Levin solve at order n.  Writes the kernel-projected
 * result into (re,im).  Returns false on sampling failure or a singular system. */
static bool levin_solve_order_mpfr(int n, double a, double b, long bits,
                                   LevinSampleMPFR amp, void* amp_ctx,
                                   LevinSampleMPFR gprime, void* gprime_ctx,
                                   LevinSampleMPFR gphase, void* gphase_ctx,
                                   LevinKernel kernel, mpfr_t re, mpfr_t im) {
    mpfr_prec_t wp = (mpfr_prec_t)bits;
    const double half = 0.5 * (b - a), mid = 0.5 * (a + b);
    bool ok = true;

    mpfr_t pi, t, xj, dtdx, sc, am, bm;
    mpfr_init2(pi, wp); mpfr_init2(t, wp); mpfr_init2(xj, wp);
    mpfr_init2(dtdx, wp); mpfr_init2(sc, wp); mpfr_init2(am, wp); mpfr_init2(bm, wp);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_set_d(dtdx, 1.0 / half, MPFR_RNDN);
    mpfr_set_d(am, a, MPFR_RNDN); mpfr_set_d(bm, b, MPFR_RNDN);

    mpfr_t* T  = malloc(sizeof(mpfr_t) * n);
    mpfr_t* dT = malloc(sizeof(mpfr_t) * n);
    ncpx* A = malloc(sizeof(ncpx) * (size_t)n * n);
    ncpx* rhs = malloc(sizeof(ncpx) * n);
    for (int k = 0; k < n; k++) { mpfr_init2(T[k], wp); mpfr_init2(dT[k], wp); }
    for (int i = 0; i < n * n; i++) ncpx_init(&A[i], wp);
    for (int j = 0; j < n; j++) ncpx_init(&rhs[j], wp);

    ncpx gpj; ncpx_init(&gpj, wp);
    for (int j = 0; j < n && ok; j++) {
        /* CGL node t_j = cos(pi j/(n-1)); x_j = mid + half t_j. */
        mpfr_mul_ui(t, pi, (unsigned long)j, MPFR_RNDN);
        mpfr_div_ui(t, t, (unsigned long)(n - 1), MPFR_RNDN);
        mpfr_cos(t, t, MPFR_RNDN);
        mpfr_mul_d(xj, t, half, MPFR_RNDN);
        mpfr_add_d(xj, xj, mid, MPFR_RNDN);
        if (!amp(amp_ctx, xj, rhs[j].re, rhs[j].im)
            || !gprime(gprime_ctx, xj, gpj.re, gpj.im)) { ok = false; break; }
        levin_cheb_row_mpfr(t, n, T, dT, wp);
        for (int k = 0; k < n; k++) {
            /* A[j,k] = dT[k]*dtdx + i g'(x_j) T[k]
               re = dT*dtdx - g'.im * T ;  im = g'.re * T */
            mpfr_mul(sc, dT[k], dtdx, MPFR_RNDN);
            mpfr_mul(A[j * n + k].im, gpj.im, T[k], MPFR_RNDN);
            mpfr_sub(A[j * n + k].re, sc, A[j * n + k].im, MPFR_RNDN);
            mpfr_mul(A[j * n + k].im, gpj.re, T[k], MPFR_RNDN);
        }
    }
    ncpx_clear(&gpj);

    if (ok) ok = levin_ncpx_solve(A, rhs, n, wp);

    if (ok) {
        /* p(b) = Σ c_k ;  p(a) = Σ (-1)^k c_k. */
        ncpx pb, pa; ncpx_init(&pb, wp); ncpx_init(&pa, wp);
        ncpx_set_ui(&pb, 0); ncpx_set_ui(&pa, 0);
        for (int k = 0; k < n; k++) {
            ncpx_add(&pb, &pb, &rhs[k]);
            if (k % 2 == 0) ncpx_add(&pa, &pa, &rhs[k]);
            else            ncpx_sub(&pa, &pa, &rhs[k]);
        }
        /* g(a), g(b); kernel factor e^{i g} = exp(i*g). */
        ncpx ga, gb, ea, eb, ta, tb, I_val;
        ncpx_init(&ga, wp); ncpx_init(&gb, wp); ncpx_init(&ea, wp);
        ncpx_init(&eb, wp); ncpx_init(&ta, wp); ncpx_init(&tb, wp); ncpx_init(&I_val, wp);
        if (gphase(gphase_ctx, am, ga.re, ga.im) && gphase(gphase_ctx, bm, gb.re, gb.im)) {
            /* i*g = (-g.im) + i (g.re) */
            mpfr_neg(ea.re, ga.im, MPFR_RNDN); mpfr_set(ea.im, ga.re, MPFR_RNDN);
            mpfr_neg(eb.re, gb.im, MPFR_RNDN); mpfr_set(eb.im, gb.re, MPFR_RNDN);
            ncpx_exp(&ea, &ea, wp); ncpx_exp(&eb, &eb, wp);
            ncpx_mul(&tb, &pb, &eb, wp);
            ncpx_mul(&ta, &pa, &ea, wp);
            ncpx_sub(&I_val, &tb, &ta);
            if (kernel == LEVIN_KERNEL_COS)      { mpfr_set(re, I_val.re, MPFR_RNDN); mpfr_set_ui(im, 0, MPFR_RNDN); }
            else if (kernel == LEVIN_KERNEL_SIN) { mpfr_set(re, I_val.im, MPFR_RNDN); mpfr_set_ui(im, 0, MPFR_RNDN); }
            else                                 { mpfr_set(re, I_val.re, MPFR_RNDN); mpfr_set(im, I_val.im, MPFR_RNDN); }
        } else ok = false;
        ncpx_clear(&ga); ncpx_clear(&gb); ncpx_clear(&ea); ncpx_clear(&eb);
        ncpx_clear(&ta); ncpx_clear(&tb); ncpx_clear(&I_val);
        ncpx_clear(&pb); ncpx_clear(&pa);
    }

    for (int k = 0; k < n; k++) { mpfr_clear(T[k]); mpfr_clear(dT[k]); }
    for (int i = 0; i < n * n; i++) ncpx_clear(&A[i]);
    for (int j = 0; j < n; j++) ncpx_clear(&rhs[j]);
    free(T); free(dT); free(A); free(rhs);
    mpfr_clear(pi); mpfr_clear(t); mpfr_clear(xj); mpfr_clear(dtdx);
    mpfr_clear(sc); mpfr_clear(am); mpfr_clear(bm);
    return ok && mpfr_number_p(re) && mpfr_number_p(im);
}

bool levin_collocation_mpfr(double a, double b, long bits,
                            LevinSampleMPFR amp,    void* amp_ctx,
                            LevinSampleMPFR gprime, void* gprime_ctx,
                            LevinSampleMPFR gphase, void* gphase_ctx,
                            LevinKernel kernel, double reltol, int n_max,
                            mpfr_t out_re, mpfr_t out_im, bool* converged) {
    mpfr_prec_t wp = (mpfr_prec_t)bits;
    if (reltol <= 0.0) reltol = 1e-10;
    if (n_max < LEVIN_N_MIN) n_max = LEVIN_N_MIN;
    *converged = false;
    if (!(b > a)) return false;

    mpfr_t pre, pim, diff, mag, tol; bool have_prev = false, have = false;
    mpfr_init2(pre, wp); mpfr_init2(pim, wp);
    mpfr_init2(diff, wp); mpfr_init2(mag, wp); mpfr_init2(tol, wp);

    for (int n = LEVIN_N_MIN; n <= n_max; n *= 2) {
        if (!levin_solve_order_mpfr(n, a, b, bits, amp, amp_ctx, gprime, gprime_ctx,
                                    gphase, gphase_ctx, kernel, out_re, out_im)) {
            if (!have) break;     /* nothing yet: not Levin-amenable */
            break;
        }
        have = true;
        if (have_prev) {
            mpfr_sub(diff, out_re, pre, MPFR_RNDN); mpfr_abs(diff, diff, MPFR_RNDN);
            mpfr_sub(mag, out_im, pim, MPFR_RNDN); mpfr_abs(mag, mag, MPFR_RNDN);
            mpfr_add(diff, diff, mag, MPFR_RNDN);        /* L1 successive change */
            mpfr_abs(mag, out_re, MPFR_RNDN);
            mpfr_abs(tol, out_im, MPFR_RNDN); mpfr_add(mag, mag, tol, MPFR_RNDN);
            mpfr_mul_d(tol, mag, reltol, MPFR_RNDN);
            if (mpfr_cmp(diff, tol) <= 0) { *converged = true; break; }
        }
        mpfr_set(pre, out_re, MPFR_RNDN); mpfr_set(pim, out_im, MPFR_RNDN);
        have_prev = true;
    }
    mpfr_clear(pre); mpfr_clear(pim); mpfr_clear(diff); mpfr_clear(mag); mpfr_clear(tol);
    return have;
}
#endif /* USE_MPFR */

LevinResult levin_collocation_machine(
    double a, double b,
    GkSampleMachine amp,    void* amp_ctx,
    GkSampleMachine gprime, void* gprime_ctx,
    GkSampleMachine gphase, void* gphase_ctx,
    LevinKernel kernel, double reltol, int n_max) {
    LevinResult r = { false, false, 0.0, INFINITY };
    if (reltol <= 0.0) reltol = 1e-10;
    if (n_max < LEVIN_N_MIN) n_max = LEVIN_N_MIN;
    if (!(b > a)) return r;   /* degenerate interval */

    double _Complex prev = 0.0;
    bool have_prev = false;
    for (int n = LEVIN_N_MIN; n <= n_max; n *= 2) {
        double _Complex In; double rcond;
        if (!levin_solve_order(n, a, b, amp, amp_ctx, gprime, gprime_ctx,
                               gphase, gphase_ctx, kernel, &In, &rcond)) {
            /* If nothing has converged yet, the problem is not Levin-amenable. */
            return r;
        }
        if (!isfinite(creal(In)) || !isfinite(cimag(In))) return r;
        r.have = true; r.val = In;
        if (have_prev) {
            double err = cabs(In - prev);
            r.err = err;
            if (err <= reltol * cabs(In)) { r.conv = true; return r; }
        }
        prev = In; have_prev = true;
    }
    return r;   /* best estimate, not converged */
}
