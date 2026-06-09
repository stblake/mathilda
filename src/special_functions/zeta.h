#ifndef ZETA_H
#define ZETA_H

#include "expr.h"

/* Zeta[s]     -- Riemann zeta function       zeta(s)   = Sum_{k>=1} k^-s
 * Zeta[s, a]  -- Hurwitz (generalized) zeta  zeta(s,a) = Sum_{k>=0} (k+a)^-s
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * Evaluation strategy (see zeta.c for detail):
 *   - Exact integer arguments reduce to closed forms: zeta(2n) to a rational
 *     multiple of Pi^(2n), zeta(-n) to a rational (both via Bernoulli numbers);
 *     zeta(0) = -1/2, zeta(1) = ComplexInfinity, odd zeta(2n+1) stays symbolic.
 *   - Exact Hurwitz at a positive integer a = m reduces to
 *     Zeta[s] - Sum_{k=1}^{m-1} k^-s.
 *   - Numeric real one-argument zeta uses MPFR's mpfr_zeta (machine or
 *     arbitrary precision).
 *   - All complex, and all two-argument (Hurwitz), numerics use one
 *     Euler-Maclaurin complex-MPFR kernel (Riemann is the a = 1 case).
 *   - Everything else stays symbolic. */
Expr* builtin_zeta(Expr* res);

void zeta_init(void);

#endif /* ZETA_H */
