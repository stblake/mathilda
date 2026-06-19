#ifndef LERCHPHI_H
#define LERCHPHI_H

#include "expr.h"

/* LerchPhi[z, s, a] -- the Lerch transcendent
 *
 *   Phi(z, s, a) = Sum_{k>=0} z^k / (k + a)^s        (|z| < 1; analytic
 *                  continuation elsewhere, branch cut z in [1, Infinity))
 *
 * For Re(a) < 0 the principal value uses the symmetric power
 * ((k + a)^2)^(-s/2), and any term with k + a = 0 is excluded. LerchPhi
 * generalizes Zeta, HurwitzZeta and PolyLog:
 *   Phi(1, s, a) = Zeta[s, a],   z Phi(z, s, 1) = PolyLog[s, z].
 *
 * Options:
 *   DoublyInfinite -> True       sum k from -Infinity to Infinity, realised as
 *                                Phi(z,s,a) + z^-1 Phi(1/z, s, 1-a).
 *   IncludeSingularTerm -> True  keep the k + a = 0 term (ComplexInfinity at a
 *                                non-positive integer a).
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * Evaluation strategy (see lerchphi.c):
 *   - Exact reductions: z = 0 -> a^-s; s = 0 -> 1/(1-z); z = 1 -> Zeta[s,a];
 *     z = -1 -> 2^-s (Zeta[s,a/2] - Zeta[s,(a+1)/2]); positive integer a ->
 *     a PolyLog form; negative integer s -> a rational function of z.
 *   - Numeric (>= 1 inexact operand): a complex-MPFR power series for |z| < 1,
 *     at machine or arbitrary precision.  |z| > 1 diverges -> stays symbolic. */
Expr* builtin_lerchphi(Expr* res);

void lerchphi_init(void);

#endif /* LERCHPHI_H */
