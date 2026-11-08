/*
 * nsum_blocks.h — summation of an integer-indexed series from its terms alone
 *
 * NSum's smooth-summand methods (Euler–Maclaurin, the tail integral) need the
 * summand off the integers.  A summand that exists only at the integers — or a
 * step function of the index such as t(Round[n]) — has no such extension, and
 * when it also converges logarithmically (1/n^2, 1/(n IntegerLength[n]^2)) the
 * partial-sum extrapolators that see only the terms (Wynn) return a confident
 * value that is wrong in the second digit.
 *
 * nsb_sum() sums Sum_{n >= n0} a(n) from integer samples, in three cases:
 *
 *   SMOOTH     the terms follow one smooth law: Levin's u-transform on the
 *              partial sums, with beta at the true index (it models c/n^p and
 *              alternating tails), error = disagreement of two lengths.
 *
 *   BLOCKED    the terms are smooth on blocks [b^(k-1)+c, b^k+c) and jump
 *              between them (anything that depends on a digit or bit length):
 *              the breakpoints are located in the head, the block structure is
 *              verified there, and the sum becomes a series over k of block
 *              sums.  Each block sum is computed by Euler–Maclaurin on the
 *              block, with the integral taken in log n by Gauss–Legendre over a
 *              local polynomial interpolant of the integer samples and the
 *              endpoint derivatives from finite differences — so a block of
 *              10^40 terms costs a few dozen samples.  The block-sum series is
 *              summed by the same routine, recursively (it may itself be
 *              blocked: 1/(n d(n) d(d(n))^2)).
 *
 *   ROUGH      anything else (1/Prime[n]^2, arithmetic functions): declined,
 *              so the caller keeps its own method.
 *
 * The error estimate is honest by construction: it is the Levin disagreement
 * plus every block's quadrature disagreement, and when a level cannot be
 * continued (index digits / time budget) the unsummed remainder is bounded
 * from the observed terms rather than dropped.
 */
#ifndef MATHILDA_NSUM_BLOCKS_H
#define MATHILDA_NSUM_BLOCKS_H

#include <stdbool.h>

#ifdef USE_MPFR
#include <gmp.h>
#include <mpfr.h>

/* The real summand at the integer index n, written to `out` at out's
 * precision.  False when it is not a real number there (the caller's routine
 * then declines). */
typedef bool (*NsbTermFn)(void* ud, const mpz_t n, mpfr_t out);

typedef struct {
    long   target_bits;   /* accuracy wanted (53 for machine precision)      */
    long   work_bits;     /* arithmetic precision (> target_bits)            */
    double budget_s;      /* wall-clock budget for the whole sum             */
} NsbParams;

/* Sum_{n >= n0} f(n).  True with *out and *err (an absolute error estimate,
 * possibly large) when the series is SMOOTH or BLOCKED; false (declined, `out`
 * untouched) when it is ROUGH, has a zero or non-real term in the head, or the
 * head itself could not be evaluated.  *levels (may be NULL) receives the
 * number of block levels used (0 for SMOOTH). */
bool nsb_sum(NsbTermFn f, void* ud, const mpz_t n0, const NsbParams* prm,
             mpfr_t out, double* err, int* levels);

#endif /* USE_MPFR */

#endif /* MATHILDA_NSUM_BLOCKS_H */
