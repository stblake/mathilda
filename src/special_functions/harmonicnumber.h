#ifndef HARMONICNUMBER_H
#define HARMONICNUMBER_H

#include "expr.h"

/* HarmonicNumber[n]     -- the n-th harmonic number  H_n     = Sum_{i=1}^n 1/i.
 * HarmonicNumber[n, r]  -- the order-r harmonic number H_n^(r) = Sum_{i=1}^n 1/i^r.
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * Evaluation strategy (see harmonicnumber.c for detail).  The function leans on
 * existing primitives (Zeta, PolyGamma, EulerGamma, BernoulliB) rather than
 * carrying its own numeric kernels:
 *   - n a non-negative integer (<= cap): expand to the explicit finite sum
 *     Sum_{i=1}^n i^-r and evaluate.  Integer r combines to an exact rational;
 *     symbolic/complex r leaves an explicit Plus (HarmonicNumber[4, r] ->
 *     1 + 2^-r + 3^-r + 4^-r).  HarmonicNumber[0, r] is 0.
 *   - n -> Infinity: HarmonicNumber[Infinity, r] = Zeta[r]
 *     (so HarmonicNumber[Infinity, 2] -> Pi^2/6; r == 1 -> ComplexInfinity).
 *   - r a non-positive integer: the Faulhaber polynomial Sum_{i=1}^n i^(-r),
 *     built from BernoulliB; symbolic n stays a polynomial in n, inexact n
 *     evaluates numerically (HarmonicNumber[z, -4] -> -z/30 + z^3/3 + ...).
 *   - otherwise, with an inexact argument and a numericizable n: reduce to the
 *     analytic form  Zeta[r] - Zeta[r, n+1]  (and, for r == 1, the equivalent
 *     EulerGamma + PolyGamma[0, n+1], avoiding the cancelling poles), then N
 *     it -- precision tracking and complex arguments come for free.
 *   - everything else stays symbolic (return NULL). */
Expr* builtin_harmonicnumber(Expr* res);

void harmonicnumber_init(void);

#endif /* HARMONICNUMBER_H */
