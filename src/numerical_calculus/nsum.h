/*
 * nsum.h — NSum[f, {i, imin, imax (, di)}, opts] and multidimensional sums
 *
 * Numerical approximation of a (finite or infinite) sum.  Mirrors Mathematica's
 * NSum: HoldAll, Block-localises the summation index, evaluates the summand
 * symbolically once per term, and accelerates the partial-sum / term sequence.
 *
 * Methods (Method -> ...):
 *   WynnEpsilon | SequenceLimit  Wynn's epsilon on the partial sums (general
 *                                fallback; excellent for alternating series).
 *   EulerMaclaurin | Integrate   explicit head terms + exp-sinh tail integral
 *                                + Bernoulli/derivative corrections (best for
 *                                monotone, slowly-converging positive series).
 *   AlternatingSigns             Cohen–Villegas–Zagier acceleration.
 *   Automatic                    chooses by sampling the summand's sign pattern.
 *
 * Options: WorkingPrecision, NSumTerms, NSumExtraTerms, WynnDegree,
 * VerifyConvergence, AccuracyGoal, PrecisionGoal (Compiled / EvaluationMonitor
 * are accepted and ignored).  Attributes: HoldAll, Protected.
 */
#ifndef MATHILDA_NSUM_H
#define MATHILDA_NSUM_H

#include "expr.h"

Expr* builtin_nsum(Expr* res);
void  nsum_init(void);

#endif /* MATHILDA_NSUM_H */
