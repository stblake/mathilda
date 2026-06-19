#ifndef HURWITZZETA_H
#define HURWITZZETA_H

#include "expr.h"

/* HurwitzZeta[s, a] -- Hurwitz zeta function
 *     zeta(s, a) = Sum_{k>=0} (k + a)^-s            (Re s > 1; elsewhere by
 *                                                    analytic continuation)
 *
 * HurwitzZeta is identical to the two-argument Zeta for Re(a) > 0, but uses a
 * different branch: each summand is the principal-branch power (k + a)^-s,
 * whereas Zeta[s, a] effectively uses ((k + a)^2)^(-s/2). The two therefore
 * disagree for non-positive real a, and HurwitzZeta keeps the singular terms
 * that Zeta drops -- so HurwitzZeta has poles at a = 0, -1, -2, ... .
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * Evaluation strategy (see hurwitzzeta.c for detail):
 *   - HurwitzZeta[s, 1] = Zeta[s] (delegates to the Riemann closed forms).
 *   - HurwitzZeta[s, 1/2] = (2^s - 1) Zeta[s].
 *   - HurwitzZeta[s, m] at a positive integer m reduces to
 *     Zeta[s] - Sum_{k=1}^{m-1} k^-s.
 *   - At a non-positive integer a: ComplexInfinity for positive integer s,
 *     the Bernoulli-polynomial value -BernoulliB[1-s, a]/(1-s) for non-positive
 *     integer s (so HurwitzZeta[0, 0] = 1/2).
 *   - Numeric (real or complex, machine or MPFR) arguments evaluate via an
 *     Euler-Maclaurin complex kernel; a numeric pole gives ComplexInfinity.
 *   - Everything else stays symbolic. */
Expr* builtin_hurwitzzeta(Expr* res);

void hurwitzzeta_init(void);

#endif /* HURWITZZETA_H */
