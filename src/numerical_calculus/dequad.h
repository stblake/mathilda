/*
 * dequad.h — double-exponential (exp-sinh) quadrature for a half line
 *
 * Evaluates  I = ∫_a^∞ f(x) dx  for an integrand that decays at +∞ and is
 * analytic on (a, ∞), via the exp-sinh transformation
 *
 *     x(t) = a + exp((π/2) sinh t),     t ∈ (-∞, ∞),
 *     I    = ∫ f(x(t)) x'(t) dt,   x'(t) = (π/2) cosh t · exp((π/2) sinh t),
 *
 * to which the trapezoidal rule converges double-exponentially.  The step h is
 * halved level by level until two successive estimates agree to the requested
 * tolerance; the tails are truncated where the transformed integrand underflows
 * relative to the running maximum.  This is the quadrature Mathematica's NSum
 * uses for the Euler–Maclaurin tail integral (Method -> DoubleExponential), and
 * is reusable groundwork for a future real-axis NIntegrate.
 *
 * The integrand is supplied through a sample callback (machine `double _Complex`
 * or an MPFR real/imag pair), mirroring quadrature.h.  A sample that returns
 * false (non-numeric / non-finite) truncates that tail rather than aborting.
 */
#ifndef MATHILDA_DEQUAD_H
#define MATHILDA_DEQUAD_H

#include <complex.h>
#include <stdbool.h>

#ifdef USE_MPFR
#  include <mpfr.h>
#endif

/* Machine sample: write f(x) into *out; return false to truncate the tail. */
typedef bool (*DeQuadSampleMachine)(void* ctx, double x, double _Complex* out);

/* Compute ∫_a^∞ f dx.  Writes *result and, via *abserr, the |last - previous|
 * level difference (an error proxy).  Returns true if the levels settled to
 * `reltol` within `max_levels` refinements, false otherwise (the best estimate
 * is still written). */
bool dequad_halfline_machine(DeQuadSampleMachine f, void* ctx, double a,
                             double reltol, int max_levels,
                             double _Complex* result, double* abserr);

#ifdef USE_MPFR
/* MPFR sample: write Re/Im of f(x) into (out_re,out_im) at the working
 * precision; return false to truncate the tail. */
typedef bool (*DeQuadSampleMPFR)(void* ctx, const mpfr_t x,
                                 mpfr_t out_re, mpfr_t out_im);

/* MPFR ∫_a^∞ f dx at `bits` precision.  (out_re,out_im) must be pre-initialised
 * by the caller.  *abserr receives the final level difference as a double. */
bool dequad_halfline_mpfr(DeQuadSampleMPFR f, void* ctx, const mpfr_t a,
                          long bits, double reltol, int max_levels,
                          mpfr_t out_re, mpfr_t out_im, double* abserr);
#endif /* USE_MPFR */

#endif /* MATHILDA_DEQUAD_H */
