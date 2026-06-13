/*
 * nlimit.h — NLimit[expr, z -> z0, opts]
 *
 * Numerically finds the limiting value of `expr` as z -> z0.  A geometric
 * sequence of sample points approaching z0 is constructed and the limit is
 * recovered by sequence acceleration / extrapolation.  Two methods:
 *
 *   Method -> EulerSum       (default) Richardson / Romberg extrapolation of
 *                            the sample sequence viewed as a function of the
 *                            geometric step.  Best for smooth (power-series)
 *                            approaches.  Depth set by Terms (default 7).
 *
 *   Method -> SequenceLimit  Wynn's epsilon algorithm (iterated Shanks
 *                            transform) on the sample sequence.  Exact in one
 *                            step for a geometric/exponential tail; controlled
 *                            by WynnDegree (default 1, needs >= 2(d+1) terms).
 *
 * Attributes: Protected (not Listable — the z -> z0 spec must not be split).
 */
#ifndef MATHILDA_NLIMIT_H
#define MATHILDA_NLIMIT_H

#include "expr.h"

Expr* builtin_nlimit(Expr* res);
void  nlimit_init(void);

#endif /* MATHILDA_NLIMIT_H */
