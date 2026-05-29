/* root_numeric.h
 * --------------
 * Numerical evaluation of Root[Function[p[#]], k] objects.
 *
 * Pipeline (in root_numeric.c):
 *   1.  Extract polynomial p from Root[Function[p_in_slot1], k].
 *   2.  Reduce to squarefree part (Yun: p / gcd(p, p')).
 *   3.  Build n*n Frobenius companion matrix at modest MPFR precision.
 *   4.  Run Hessenberg + Francis QR via `eigen_all_eigenvalues_real_mpfr`
 *       to obtain all n approximate roots.
 *   5.  Use Sturm sign-variation to certify the real-root count; escalate
 *       precision if the QR classification disagrees.
 *   6.  Canonically sort the roots (real-first ascending; complex by
 *       Re ascending, |Im| ascending, smaller-Im first within a
 *       conjugate pair).
 *   7.  Pick the k-th candidate (1-indexed) and refine via fused
 *       Horner-Newton at full requested precision (plus guard bits).
 *   8.  Verify the refined root has not jumped basins.
 *   9.  Return EXPR_MPFR (real root) or Complex[EXPR_MPFR, EXPR_MPFR]
 *       (complex root).
 *
 * The shared canonical comparator is exposed for `solvepoly.c` so the
 * `k` emitted by `Solve` agrees with the `k` interpreted by `N`.
 */

#ifndef ROOT_NUMERIC_H
#define ROOT_NUMERIC_H

#include "expr.h"
#include "numeric.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Numerically evaluate `root_expr` (must be `Root[Function[..], k]`) at
 * the precision encoded in `spec`. Returns a freshly allocated Expr
 * (EXPR_MPFR or Complex[EXPR_MPFR, EXPR_MPFR] under USE_MPFR; EXPR_REAL
 * or Complex[EXPR_REAL, EXPR_REAL] under machine precision), or NULL if
 * the input cannot be numericalized (non-integer coefficients, out-of-
 * range k, QR non-convergence, etc.). On NULL the caller should pass
 * the original Root expression through unchanged.
 *
 * Diagnostics are printed to stderr in the `Root::<tag>` style used by
 * FindRoot. */
Expr* root_numericalize(const Expr* root_expr, NumericSpec spec);

#ifdef USE_MPFR
/* Canonical Mathematica root ordering, comparing two complex MPFR
 * candidates carrying a real/complex classification.
 *
 *   real-first: real (true) < complex (false).
 *   reals:      ascending by re.
 *   complex:    ascending by re; tie-break by |im| ascending; final tie-
 *               break by sign(im) (negative first, so a conjugate pair
 *               (a - bi, a + bi) sorts as (a - bi, a + bi)).
 *
 * Used by both root_numericalize and solvepoly.c. Returns < 0, 0, > 0. */
int root_canonical_cmp_mpfr(const mpfr_t a_re, const mpfr_t a_im, int a_real,
                            const mpfr_t b_re, const mpfr_t b_im, int b_real);
#endif

#ifdef __cplusplus
}
#endif

#endif /* ROOT_NUMERIC_H */
