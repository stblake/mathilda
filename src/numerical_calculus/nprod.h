/*
 * nprod.h — NProduct[f, {i, imin, imax (, di)}, opts] and multidimensional
 * products.
 *
 * Numerical approximation of a (finite or infinite) product.  Mirrors
 * Mathematica's NProduct: HoldAll, Block-localises the product index, and — per
 * Keiper 1992 — is evaluated as Exp[NSum[Log[f], ...]], reusing the full NSum
 * engine (Euler-Maclaurin / Wynn epsilon / Cohen-Villegas-Zagier, MPFR,
 * divergence detection, multidimensional recursion).
 *
 * Options: Method, WorkingPrecision, NProductFactors, NProductExtraFactors,
 * WynnDegree, VerifyConvergence, AccuracyGoal, PrecisionGoal (EvaluationMonitor
 * accepted and ignored).  Attributes: HoldAll, Protected.
 */
#ifndef MATHILDA_NPROD_H
#define MATHILDA_NPROD_H

#include "expr.h"

Expr* builtin_nprod(Expr* res);
void  nprod_init(void);

#endif /* MATHILDA_NPROD_H */
