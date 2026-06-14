/*
 * nint.h — NIntegrate[f, {x, xmin, xmax (, x1, ...)}, ...] numerical integration
 *
 * Numerical approximation of a definite integral.  Mirrors Mathematica's
 * NIntegrate: HoldAll, Block-localises the integration variable(s), evaluates
 * the integrand symbolically at numeric sample points, and drives a suite of
 * quadrature engines selected by Method / by analysis of the region:
 *
 *   GlobalAdaptive / GaussKronrodRule  adaptive 7-15 Gauss-Kronrod (finite, smooth)
 *   DoubleExponential (tanh-sinh / sinh-sinh / exp-sinh)  endpoint singularities,
 *                                       infinite ranges, arbitrary precision
 *   LevinRule                          oscillatory integrands (zeros + Wynn / Levin)
 *   MonteCarlo / QuasiMonteCarlo / AdaptiveMonteCarlo  high dimensions, regions
 *   PrincipalValue                     Cauchy principal values about Exclusions
 *
 * Multidimensional integrals are evaluated by iterated 1D quadrature (the
 * integrand of the outer variable is an inner NIntegrate over the rest), so an
 * inner bound may depend on an outer variable.  Complex xmin/xmax (or extra
 * intermediate nodes) define a piecewise-linear contour in the complex plane.
 *
 * Options: WorkingPrecision (default MachinePrecision), PrecisionGoal,
 * AccuracyGoal, MaxRecursion, MinRecursion, MaxPoints, Method, Exclusions
 * (EvaluationMonitor accepted and ignored).  Attributes: HoldAll, Protected.
 *
 * Memory: receives `res` owned by the evaluator; returns a fresh Expr* on
 * success or NULL (unevaluated).  Never frees `res`.  Every temporary variable
 * binding is removed on all return paths.
 */
#ifndef MATHILDA_NINT_H
#define MATHILDA_NINT_H

#include "expr.h"

Expr* builtin_nintegrate(Expr* res);
void  nintegrate_init(void);

#endif /* MATHILDA_NINT_H */
