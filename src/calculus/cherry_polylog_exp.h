/* cherry_polylog_exp.h — exponential-tower polylogarithm-ladder integration.
 *
 * The general-weight, algebraic-root generalisation of rt_cherry_dilog_exp
 * (cherry_dilog_exp.c): integrands P(x)/Q(E^(c x)) that are rational in the
 * exponential kernel theta = E^(c x) with polynomial-in-x coefficients.  Partial
 * fractioning over the (rational OR algebraic) roots rho of Q reduces the problem
 * to simple poles x^n/(theta - rho), each closed by the exact Cherry polylogarithm
 * ladder
 *
 *   INT x^n/(theta - rho) dx
 *     = Sum_{k=0}^{n} -(1/rho) (n!/(n-k)!) / c^(k+1) x^(n-k) PolyLog[k+1, rho/theta]
 *
 * (PolyLog[1, z] = -Log[1-z]), so an x^n numerator yields polylogarithms up to
 * weight n+1.  Covers e.g.
 *   INT x^2/(E^x-1) dx        = x^2 Log[1-E^-x] - 2 x PolyLog[2,E^-x]
 *                               - 2 PolyLog[3,E^-x]
 *   INT x^4/(E^(5x)-1) dx     = ... - 24/3125 PolyLog[5, E^-5x]
 *   INT x/(E^(2x)+E^x-1) dx   (algebraic roots (-1+-Sqrt[5])/2)
 *
 * Returns a fresh, PowerExpand-diff-back-verified antiderivative or NULL.
 * Defined in cherry_polylog_exp.c.
 */

#ifndef MATHILDA_CHERRY_POLYLOG_EXP_H
#define MATHILDA_CHERRY_POLYLOG_EXP_H

#include "expr.h"

Expr* rt_cherry_polylog_exp(Expr* f, Expr* x);

/* Integrate`Cherry`PolyLogExp[f, x] — direct debuggable surface (bypasses the
 * cascade).  Registered by cherry_builtins_init (cherry_driver.c). */
Expr* builtin_cherry_polylog_exp(Expr* res);

#endif /* MATHILDA_CHERRY_POLYLOG_EXP_H */
