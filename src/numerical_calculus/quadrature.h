/*
 * quadrature.h — periodic-trapezoidal contour integration core.
 *
 * A reusable, expression-agnostic numerical engine for evaluating the
 * Cauchy contour integral
 *
 *     Res(f, z0) = (1/2pi i) oint_{|z-z0|=r} f(z) dz
 *                = (r/N) sum_{k=0}^{N-1} f(z0 + r e^{i theta_k}) e^{i theta_k}
 *
 * by the periodic trapezoidal rule (theta_k = 2 pi k / N). For a function
 * analytic on the integration circle the rule converges geometrically
 * (spectral accuracy — Trefethen & Weideman, SIAM Review 2014), so a small
 * number of N-doublings reaches full precision.
 *
 * The engine knows nothing about Mathilda expressions: the caller supplies
 * a sampler callback that returns f(z) as a raw complex number (machine
 * `double _Complex`, or an MPFR real/imag pair). This keeps the core small,
 * independently testable, and reusable by a future NIntegrate.
 *
 * Features beyond a plain fixed-N rule:
 *   - Adaptive N-doubling with sample reuse, stopping at a precision goal.
 *   - Aitken/Shanks (dea3-style) extrapolation of the doubling sequence,
 *     squeezing extra digits and supplying an error estimate.
 *   - Adaptive optimal radius (Fornberg/Bornemann style) when the caller
 *     requests it (radius <= 0): the true residue is radius-independent, so
 *     the search simply favours radii that converge fastest.
 *   - Branch-cut / non-convergence detection from the decay rate of the
 *     doubling errors (a non-geometric tail flags a non-analytic contour).
 */
#ifndef MATHILDA_QUADRATURE_H
#define MATHILDA_QUADRATURE_H

#include <complex.h>
#include <stdbool.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Outcome of a contour-quadrature run. */
typedef enum {
    QD_OK = 0,      /* converged to the requested tolerance              */
    QD_NOCONV,      /* converging but too slow (ncvi-style)              */
    QD_BRANCHCUT,   /* errors did not decay geometrically (non-analytic) */
    QD_NONNUMERIC   /* sampler could not return a finite number on contour */
} QdStatus;

/* Machine-precision sampler. Given z = z0 + r e^{i theta}, write f(z) into
 * *out and return true; return false (or a non-finite value) to signal that
 * f(z) is not a usable number there. `ctx` is the caller's opaque cookie. */
typedef bool (*QdSampleMachine)(void* ctx, double _Complex z,
                                double _Complex* out);

typedef struct {
    double _Complex value;   /* estimated residue                    */
    double          abs_err; /* estimated absolute error             */
    QdStatus        status;
    int             n_used;  /* sample points at termination         */
    double          radius;  /* radius actually used                 */
} QdResultMachine;

/* Machine-precision contour residue.
 *   radius > 0  : use this fixed radius.
 *   radius <= 0 : Automatic — search for a fast-converging radius.
 *   prec_goal_digits : target accuracy, in decimal digits.
 *   max_recursion    : maximum number of N-doublings (>= 0). */
QdResultMachine qd_contour_residue_machine(QdSampleMachine f, void* ctx,
                                           double _Complex z0,
                                           double radius,
                                           double prec_goal_digits,
                                           int max_recursion);

#ifdef USE_MPFR
/* MPFR sampler. (z_re, z_im) is the sample point; write f(z) into the
 * already-init2'd (out_re, out_im). Return false on failure. */
typedef bool (*QdSampleMpfr)(void* ctx,
                             const mpfr_t z_re, const mpfr_t z_im,
                             mpfr_t out_re, mpfr_t out_im);

typedef struct {
    QdStatus status;
    double   abs_err;  /* error gauge (a double is enough as a magnitude) */
    int      n_used;
    double   radius;
} QdResultMpfr;

/* MPFR contour residue. `out_re`/`out_im` must be init2'd at `bits` by the
 * caller and receive the result. Semantics of radius / prec_goal_digits /
 * max_recursion match the machine version. */
QdResultMpfr qd_contour_residue_mpfr(QdSampleMpfr f, void* ctx,
                                     const mpfr_t z0_re, const mpfr_t z0_im,
                                     double radius, long bits,
                                     double prec_goal_digits,
                                     int max_recursion,
                                     mpfr_t out_re, mpfr_t out_im);
#endif /* USE_MPFR */

#endif /* MATHILDA_QUADRATURE_H */
