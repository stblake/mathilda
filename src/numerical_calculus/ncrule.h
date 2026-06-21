/*
 * ncrule.h — equally-spaced composite quadrature rules (Riemann / Trapezoidal /
 *            Newton–Cotes), machine and arbitrary (MPFR) precision
 *
 * These are the "fixed sampling" Newton–Cotes family of NIntegrate strategies,
 * the counterparts to the adaptive Gauss–Kronrod and double-exponential kernels:
 *
 *   RiemannRule       composite left / right / midpoint rectangle rule (first
 *                     order; midpoint second order).  No extrapolation — this is
 *                     the raw Riemann sum, refined by panel doubling.
 *
 *   TrapezoidalRule   composite trapezoidal rule.  With Romberg extrapolation
 *                     (the default) successive halvings feed a Richardson table
 *                     (classic Romberg integration); without it, the plain
 *                     piecewise-linear estimate is refined by doubling.
 *
 *   NewtonCotesRule   composite closed Newton–Cotes rule of a chosen order
 *                     (Points = 2 trapezoid, 3 Simpson, 4 Simpson-3/8, 5 Boole),
 *                     refined by doubling with Romberg (Richardson) acceleration.
 *
 * Unlike Gauss–Kronrod these rules sample the interval endpoints, matching the
 * Wolfram-Language semantics of the corresponding methods; an integrand that is
 * non-numeric at an endpoint therefore aborts the rule (no estimate) rather than
 * being silently skipped.
 *
 * The integrand is supplied through the same sample-callback types used by the
 * other kernels: GkSampleMachine (machine `double _Complex`) and, for the MPFR
 * variant, DeQuadSampleMPFR (an MPFR real/imag pair).  Each refinement compares
 * the current best estimate with the previous one; *abserr receives the final
 * level difference and the return value reports whether it settled to the
 * requested tolerance within the level / evaluation budget.
 */
#ifndef MATHILDA_NCRULE_H
#define MATHILDA_NCRULE_H

#include <complex.h>
#include <stdbool.h>

#include "gkadapt.h"   /* GkSampleMachine */
#include "dequad.h"    /* DeQuadSampleMPFR (under USE_MPFR) */

#ifdef USE_MPFR
#  include <mpfr.h>
#endif

typedef enum {
    NCR_RIEMANN_LEFT = 0,   /* left rectangle rule                              */
    NCR_RIEMANN_RIGHT,      /* right rectangle rule                             */
    NCR_RIEMANN_MIDPOINT,   /* midpoint rule                                    */
    NCR_TRAPEZOIDAL,        /* trapezoidal rule (optionally Romberg)            */
    NCR_NEWTONCOTES         /* closed Newton–Cotes of order nc_points-1         */
} NcrRuleKind;

/* Compute ∫_a^b f dx with the requested fixed rule (machine precision).
 *   nc_points   points per elementary Newton–Cotes panel (2..5); ignored for
 *               Riemann / Trapezoidal.
 *   romberg     enable Richardson (Romberg) extrapolation over the halvings
 *               (honoured for Trapezoidal / Newton–Cotes; Riemann is never
 *               extrapolated).
 *   abstol,reltol  termination: |Iⱼ − Iⱼ₋₁| <= max(abstol, reltol·|Iⱼ|).
 *   max_levels  cap on panel-doubling refinements.
 *   max_eval    cap on integrand evaluations (<=0 => unlimited).
 * Writes *result and, via *abserr, the final level difference.  Returns true if
 * it converged to the tolerance; false otherwise (best estimate still written),
 * and false with *result untouched if the integrand was non-numeric. */
bool ncr_integrate_machine(GkSampleMachine f, void* ctx, double a, double b,
                           NcrRuleKind kind, int nc_points, bool romberg,
                           double reltol, double abstol,
                           int max_levels, long max_eval,
                           double _Complex* result, double* abserr);

#ifdef USE_MPFR
/* Arbitrary-precision counterpart at `bits` working precision.  (out_re,out_im)
 * must be pre-initialised by the caller.  *abserr receives the final level
 * difference as a double. */
bool ncr_integrate_mpfr(DeQuadSampleMPFR f, void* ctx,
                        const mpfr_t a, const mpfr_t b,
                        NcrRuleKind kind, int nc_points, bool romberg, long bits,
                        double reltol, double abstol,
                        int max_levels, long max_eval,
                        mpfr_t out_re, mpfr_t out_im, double* abserr);
#endif /* USE_MPFR */

#endif /* MATHILDA_NCRULE_H */
