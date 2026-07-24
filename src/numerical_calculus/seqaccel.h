/*
 * seqaccel.h — shared sequence-acceleration kernels
 *
 * Pure-numeric extrapolation of a finite sequence S[0..terms-1] to its limit.
 * Two classic accelerators, each in a machine (`double _Complex`) and an MPFR
 * variant:
 *
 *   Richardson / Romberg ("EulerSum"): the all-powers tableau with denominator
 *     2^j - 1,
 *         T(i,0) = S_i,  T(i,j) = T(i,j-1) + (T(i,j-1) - T(i-1,j-1))/(2^j - 1),
 *     result = T(terms-1, terms-1).  Best for a geometric/power-series step.
 *
 *   Wynn's epsilon ("SequenceLimit"): the iterated Shanks transform.  The
 *     degree-d estimate lives in column 2d (stored shifted by one); the entry
 *     that best agrees with its neighbour is returned to avoid amplifying the
 *     tiny Shanks denominators near the bottom corner.
 *
 *   Levin's transformation ("Levin", Levin 1972): a nonlinear transform driven
 *     by explicit remainder estimates omega_i formed from the increments
 *     a_i = S_i - S_{i-1}.  Three variants select omega (see SeqaccelLevinVariant):
 *       t: omega_i = a_i        u: omega_i = (beta+i) a_i
 *       v: omega_i = a_i a_{i+1} / (a_i - a_{i+1})
 *     The degree-k estimate is the direct ratio of sums
 *         L = [ sum_j (-1)^j C(k,j) ((beta+n+j)/(beta+n+k))^{k-1} S_{n+j}/omega_{n+j} ]
 *           / [ sum_j (-1)^j C(k,j) ((beta+n+j)/(beta+n+k))^{k-1}      1/omega_{n+j} ]
 *     with base index n=1 (a_i needs i>=1); the power ratio is bounded <=1, so the
 *     direct form is numerically clean at the small degrees these callers use.
 *     Excels on logarithmically/algebraically convergent sequences where the
 *     power-tail Richardson tableau is weak.
 *
 * These are deliberately Expr-agnostic and apply no acceptance/noise policy:
 * each writes the extrapolated value and, via *step, the magnitude of the last
 * extrapolation step so the caller can run its own convergence gate.  Shared by
 * NLimit (sample sequences approaching a limit point) and NSum (partial-sum
 * sequences of an infinite series).
 */
#ifndef MATHILDA_SEQACCEL_H
#define MATHILDA_SEQACCEL_H

#include <complex.h>
#include <stdbool.h>

#ifdef USE_MPFR
#  include <mpfr.h>
#endif

/* Richardson/Romberg.  Writes *result and *step (|last diagonal - previous|).
 * Returns false only on allocation failure (terms >= 1 required). */
bool seqaccel_richardson_machine(const double _Complex* S, int terms,
                                 double _Complex* result, double* step);

/* Wynn's epsilon.  `degree` is clamped to [1, (terms-1)/2].  Writes *result
 * and *step (best neighbour disagreement in the result column).  Returns false
 * if terms < 3 (no Shanks step possible) or on allocation failure. */
bool seqaccel_wynn_machine(const double _Complex* S, int terms, int degree,
                           double _Complex* result, double* step);

/* Levin transformation remainder-estimate variant. */
typedef enum {
    SEQACCEL_LEVIN_U = 0,   /* omega_i = (beta+i) (S_i - S_{i-1})  */
    SEQACCEL_LEVIN_T,       /* omega_i = S_i - S_{i-1}             */
    SEQACCEL_LEVIN_V        /* omega_i = a_i a_{i+1}/(a_i - a_{i+1}) */
} SeqaccelLevinVariant;

/* Levin's transformation.  `variant` is a SeqaccelLevinVariant; `beta` is the
 * shift parameter (1.0 is standard).  Writes *result and *step (|L_k - L_{k-1}|,
 * the disagreement between the two deepest degrees).  Returns false if terms is
 * too small (u/t need >= 3, v needs >= 4), on a degenerate/zero remainder
 * estimate, or on allocation failure. */
bool seqaccel_levin_machine(const double _Complex* S, int terms, int variant,
                            double beta, double _Complex* result, double* step);

#ifdef USE_MPFR
/* MPFR variants.  (out_re, out_im) must be pre-initialised by the caller; the
 * extrapolated value is written there.  *step receives the L1 step magnitude as
 * a double (for the caller's gate) and *finite whether the result is a finite
 * number.  Return false only when the tableau cannot be built (terms too small
 * / allocation failure), in which case (out_re,out_im) are left untouched. */
bool seqaccel_richardson_mpfr(const mpfr_t* Sr, const mpfr_t* Si, int terms,
                              long bits, mpfr_t out_re, mpfr_t out_im,
                              double* step, bool* finite);

bool seqaccel_wynn_mpfr(const mpfr_t* Sr, const mpfr_t* Si, int terms,
                        int degree, long bits, mpfr_t out_re, mpfr_t out_im,
                        double* step, bool* finite);

bool seqaccel_levin_mpfr(const mpfr_t* Sr, const mpfr_t* Si, int terms,
                         int variant, double beta, long bits,
                         mpfr_t out_re, mpfr_t out_im, double* step, bool* finite);
#endif /* USE_MPFR */

#endif /* MATHILDA_SEQACCEL_H */
