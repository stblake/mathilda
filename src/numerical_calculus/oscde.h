/*
 * oscde.h — Ooura–Mori double-exponential quadrature for Fourier-type integrals
 *
 * Evaluates a semi-infinite oscillatory integral
 *
 *     I = ∫_0^∞ g(x) sin(ω x) dx     or     ∫_0^∞ g(x) cos(ω x) dx
 *
 * with a slowly-decaying amplitude g, using the Ooura–Mori (1999) double
 * exponential transformation
 *
 *     x(t) = (M/ω) φ(t),   φ(t) = t / (1 - exp(-2t - α(1-e^{-t}) - β(e^t-1))),
 *     M = π/h,   β = 1/4,   α = β / sqrt(1 + M ln(1+M)/(2π)),
 *
 * to which the (offset) trapezoidal rule is applied.  As t → +∞, φ(t) → t, so
 * the abscissae ω x(kh) → kπ approach the zeros of the oscillation double
 * exponentially — the oscillatory tail is annihilated with a handful of nodes
 * per decade rather than one panel per half-period.  As t → -∞, φ'(t) → 0
 * double exponentially, taming the x = 0 end.  This is the method that makes a
 * high-precision Fourier-type integral (∫_0^∞ x Sin[x]/(x²+4) dx, Bessel tails,
 * Sin[x]/x, …) cheap where integrate-between-the-zeros + sequence acceleration
 * is not.  Reference: T. Ooura and M. Mori, "A robust double exponential
 * formula for Fourier-type integrals", J. Comput. Appl. Math. 112 (1999).
 *
 * The whole integrand (amplitude times the trigonometric factor) is supplied
 * through the same DeQuadSampleMPFR callback used by the other DE kernels.
 */
#ifndef MATHILDA_OSCDE_H
#define MATHILDA_OSCDE_H

#include <stdbool.h>

#ifdef USE_MPFR
#include <mpfr.h>
#include "dequad.h"   /* DeQuadSampleMPFR */

/* Kind of oscillation, selecting the node grid (integer for sine so a node
 * lands on every zero kπ/ω; half-integer for cosine, zeros at (k+½)π/ω). */
typedef enum { OSCDE_SIN, OSCDE_COS } OscDeKind;

/* Compute I = ∫_0^∞ f(x) dx where f oscillates as {sin|cos}(ω x) (ω > 0) at
 * `bits` precision.  (out_re,out_im) must be pre-initialised by the caller.
 * *abserr receives the final level difference as a double.  Returns true if the
 * step-halving settled to `reltol`, false otherwise (the best estimate is still
 * written, so a caller may fall back to another method on false). */
bool oscde_fourier_mpfr(DeQuadSampleMPFR f, void* ctx, double omega,
                        OscDeKind kind, long bits, double reltol,
                        int max_levels, mpfr_t out_re, mpfr_t out_im,
                        double* abserr);
#endif /* USE_MPFR */

#endif /* MATHILDA_OSCDE_H */
