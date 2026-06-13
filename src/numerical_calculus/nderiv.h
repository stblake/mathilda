/*
 * nderiv.h — ND[expr, x, x0] / ND[expr, {x, n}, x0, opts]
 *
 * Numerical approximation to the (n-th) derivative of `expr` with respect to
 * `x` at the point x = x0. Two methods:
 *
 *   Method -> EulerSum   (default) Richardson extrapolation of forward,
 *                        direction-Scale finite differences. Works for
 *                        non-analytic expr; needs integer order n >= 1.
 *   Method -> NIntegrate Cauchy integral formula, computed by reusing
 *                        NResidue: f^(n)(x0) = Gamma(n+1) * NResidue[
 *                        expr/(x-x0)^(n+1), {x, x0}, Radius -> Scale].
 *                        Requires expr analytic near x0; supports
 *                        fractional / complex order.
 *
 * Attributes: Protected (list-threading over arg 1 is handled manually so the
 * {x, n} spec is never split — see nderiv.c).
 */
#ifndef MATHILDA_NDERIV_H
#define MATHILDA_NDERIV_H

#include "expr.h"

Expr* builtin_nd(Expr* res);
void  nd_init(void);

#endif /* MATHILDA_NDERIV_H */
