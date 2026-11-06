/*
 * denint.h — double-exponential quadrature on finite and doubly-infinite ranges
 *
 * Companion to dequad.h (which handles the half line ∫_a^∞).  Two transforms:
 *
 *   tanh-sinh, finite [a,b]:
 *     x(t) = (a+b)/2 + (b-a)/2 · tanh((π/2) sinh t),
 *     the standard double-exponential rule.  The node density collapses
 *     double-exponentially toward both endpoints, so integrable algebraic or
 *     logarithmic endpoint singularities are handled without special casing —
 *     the abscissae are formed from the accurate endpoint offset (not by
 *     evaluating tanh and subtracting) so a sample never lands exactly on a
 *     singular endpoint until its weight has already underflowed.
 *
 *   sinh-sinh, whole line (-∞, ∞):
 *     x(t) = sinh((π/2) sinh t),  for integrands decaying at both ends.
 *
 * As in dequad.h the integrand is supplied through a sample callback (the same
 * DeQuadSampleMachine / DeQuadSampleMPFR types).  A sample that returns false
 * truncates that tail rather than aborting.  Each refinement halves the step
 * and compares to the previous estimate; *abserr receives the final level
 * difference and the return value is whether it settled to `reltol`.
 */
#ifndef MATHILDA_DENINT_H
#define MATHILDA_DENINT_H

#include <complex.h>
#include <stdbool.h>

#include "dequad.h"   /* DeQuadSampleMachine / DeQuadSampleMPFR */

#ifdef USE_MPFR
#  include <mpfr.h>
#endif

/* tanh-sinh on a finite interval [a,b] (machine). */
bool denint_tanhsinh_machine(DeQuadSampleMachine f, void* ctx, double a, double b,
                             double reltol, int max_levels,
                             double _Complex* result, double* abserr);

/* sinh-sinh on the whole real line (machine). */
bool denint_sinhsinh_machine(DeQuadSampleMachine f, void* ctx,
                             double reltol, int max_levels,
                             double _Complex* result, double* abserr);

#ifdef USE_MPFR
/* tanh-sinh on [a,b] at `bits` precision; (out_re,out_im) pre-initialised. */
bool denint_tanhsinh_mpfr(DeQuadSampleMPFR f, void* ctx,
                          const mpfr_t a, const mpfr_t b,
                          long bits, double reltol, int max_levels,
                          mpfr_t out_re, mpfr_t out_im, double* abserr);

/* The same with an absolute tolerance: the levels also count as settled once
 * their difference is <= |abstol| (NULL: relative test only).  Used for a
 * piece of a split interval whose share of the whole is small. */
bool denint_tanhsinh_mpfr_abs(DeQuadSampleMPFR f, void* ctx,
                              const mpfr_t a, const mpfr_t b,
                              long bits, double reltol, mpfr_srcptr abstol,
                              int max_levels,
                              mpfr_t out_re, mpfr_t out_im, double* abserr);

/* Clenshaw-Curtis on [a,b] at `bits` precision, N = 8, 16, ..., 256 nested
 * Chebyshev points, for an integrand analytic on the interval's own scale
 * (geometric convergence; tanh-sinh suits endpoint singularities).  `cache`
 * holds the node/weight tables, built lazily and shared by the calls that
 * pass it (ccquad_cache_new(bits) / ccquad_cache_free).  Returns true and
 * writes the result only when it converged (else nothing is written: the
 * caller falls back to tanh-sinh). */
#define CC_N0 8
#define CC_MAX_LEVEL 5
typedef struct CcCache CcCache;
CcCache* ccquad_cache_new(long bits);
void     ccquad_cache_free(CcCache* c);
bool ccquad_mpfr(DeQuadSampleMPFR f, void* ctx, const mpfr_t a, const mpfr_t b,
                 long bits, double reltol, mpfr_srcptr abstol, CcCache* cache,
                 mpfr_t out_re, mpfr_t out_im, double* abserr);
/* As ccquad_mpfr, stopping at N = CC_N0 * 2^max_level (clamped to CC_MAX_LEVEL):
 * a bounded probe, which spends at most N + 1 samples before it declines, and
 * declines early (from N = 32 on) once a doubling gains under 8 bits. */
bool ccquad_mpfr_upto(DeQuadSampleMPFR f, void* ctx, const mpfr_t a, const mpfr_t b,
                      long bits, double reltol, mpfr_srcptr abstol, CcCache* cache,
                      int max_level, mpfr_t out_re, mpfr_t out_im, double* abserr);

/* sinh-sinh on the whole line at `bits` precision. */
bool denint_sinhsinh_mpfr(DeQuadSampleMPFR f, void* ctx,
                          long bits, double reltol, int max_levels,
                          mpfr_t out_re, mpfr_t out_im, double* abserr);

/* Recentred double-exponential rule, for an integrand whose mass sits in a
 * feature far narrower than, or far from the centre of, what the standard
 * rule resolves (a narrow peak exp(-x^2 - b/x^2), x^n e^-x for large n, a
 * Lorentzian 1/((x - x0)^2 + eps^2), a near-pole 1/(eps^2 + x^2) at an
 * endpoint).
 *
 * Both standard rules are a trapezoidal sum in t after the map
 * y = W sinh t, where y is a log variable of the range: y = log(x - a) on the
 * half line (exp-sinh, W = pi/2) and y = log((x - a)/(b - x)) on [a, b]
 * (tanh-sinh, W = pi).  They therefore resolve y on a scale of about
 * sqrt(y^2 + W^2) h: a feature of width w in y at y0 needs a step
 * h ~ w / |y0|, and the node count grows as 1/h -- 10^4 samples for the
 * peak of exp(-x^2 - 10^6/x^2), whose width is 1% of its abscissa.
 * Centring the same map on the feature, y = y0 + w sinh v, gives it a width
 * of about 1 in v, and the trapezoidal sum in v converges double-
 * exponentially again, now with h ~ 1/8 (the rest of the range still decays
 * double-exponentially in v, as before in t).
 *
 * denint_recentre_locate_mpfr finds (y0, w): the maximum of the log mass
 * log|f(x(y)) x'(y)| on a coarse grid in y, refined by golden-section search
 * until the bracket lies within the peak's top half-bit, and w from the
 * curvature there.  It returns true only when centring pays: when the
 * feature is narrower than the standard rule's own resolution at y0 by a
 * factor >= 4.  An oscillatory or multi-modal integrand (no dominant
 * maximum, or one at the edge of the grid) is declined.  The cost is the
 * grid (23-33 samples) plus ~log1.6(4/w) golden-section samples.
 *
 * denint_recentre_mpfr then integrates with the nested trapezoidal driver of
 * the tanh-sinh rule (same tail guard and convergence test).  b is ignored
 * for the half line; (out_re, out_im) are pre-initialised. */
typedef enum { DEN_RC_HALFLINE = 0, DEN_RC_FINITE = 1 } DenRecentreKind;
bool denint_recentre_locate_mpfr(DeQuadSampleMPFR f, void* ctx, DenRecentreKind kind,
                                 const mpfr_t a, mpfr_srcptr b, long bits,
                                 double* y0, double* w);
bool denint_recentre_mpfr(DeQuadSampleMPFR f, void* ctx, DenRecentreKind kind,
                          const mpfr_t a, mpfr_srcptr b, long bits,
                          double y0, double w, double reltol, int max_levels,
                          mpfr_t out_re, mpfr_t out_im, double* abserr);
#endif /* USE_MPFR */

#endif /* MATHILDA_DENINT_H */
