#ifndef POLYLOG_H
#define POLYLOG_H

#include "expr.h"

/* PolyLog[n, z]     -- polylogarithm  Li_n(z) = Sum_{k>=1} z^k / k^n
 * PolyLog[n, p, z]  -- Nielsen generalized polylogarithm S_{n,p}(z)
 *                      (accepted; left symbolic -- no closed-form engine)
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * Evaluation strategy (see polylog.c for detail):
 *   - Exact special cases reduce in closed form:
 *       PolyLog[n, 0]  = 0
 *       PolyLog[1, z]  = -Log[1-z]
 *       PolyLog[0, z]  = z/(1-z)
 *       PolyLog[-m, z] = Eulerian-number rational function (m >= 1)
 *       PolyLog[n, 1]  = Zeta[n]                 (integer n >= 2)
 *       PolyLog[n, -1] = (2^(1-n) - 1) Zeta[n]   (integer n >= 2)
 *       PolyLog[2, 1/2], PolyLog[3, 1/2] -- famous closed forms
 *   - Numeric (machine or arbitrary-precision MPFR, real or complex):
 *       |z| <= 1/2          : direct power series
 *       1/2 < |z|, |ln z|<2pi: Jonquiere / zeta expansion
 *   - Everything else stays symbolic. */
Expr* builtin_polylog(Expr* res);

void polylog_init(void);

#endif /* POLYLOG_H */
