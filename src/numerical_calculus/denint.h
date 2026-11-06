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
#endif /* USE_MPFR */

#endif /* MATHILDA_DENINT_H */
