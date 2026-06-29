/* integrate_chebychev.h
 *
 * Integration of Chebychev binomial differentials -- integrands of the form
 *
 *     x^p (a x^r + b)^q ,   p, q, r rational,  a, b free of x .
 *
 * Chebychev (1853) proved this is elementary if and only if at least one of
 *
 *     q ,   (p + 1)/r ,   q + (p + 1)/r
 *
 * is an integer, each case admitting a substitution that turns the integrand
 * into a *rational function* of a new variable:
 *
 *   Type I   (q in Z):           x = u^N, N = LCM(den p, den r).
 *   Type II  ((p+1)/r in Z):     u^s = a x^r + b, s = den q.
 *   Type III (q+(p+1)/r in Z):   u = x^r, then t^s = (a u + b)/u, s = den q.
 *
 * The reduced rational integral is handled by the BronsteinRational stage and
 * the result is back-substituted.  Recognition is a single structural scan, so
 * the method is cheap enough to run ahead of the Eliminate/Solve search in
 * DerivativeDivides.  Non-elementary integrands (none of the three quantities
 * integral) return NULL, so the Automatic cascade falls through to the
 * remaining methods.
 *
 * Examples:
 *     Integrate[x^(1/3) (a Sqrt[x] + b)^3, x]           (Type I)
 *     Integrate[(a x^2 + b)^(1/3)/x^3 ... ]             (Type II / III)
 *     Integrate[x^(-5/6) Sqrt[x^(1/3) - 1], x]          (Type III, branch-fixed)
 */
#ifndef INTEGRATE_CHEBYCHEV_H
#define INTEGRATE_CHEBYCHEV_H

#include "expr.h"

/* Direct entry point for the Integrate dispatcher's Automatic cascade.
 * Recognises x^p (a x^r + b)^q, performs the appropriate Chebychev
 * substitution, recurses into the full integrator on the rationalised
 * integrand, and back-substitutes.  Returns a freshly owned antiderivative of
 * `f` with respect to symbol `x`, or NULL when `f` is not a Chebychev binomial,
 * the integral is non-elementary, or the reduced rational integral does not
 * close.  Does NOT take ownership of `f` or `x`. */
Expr* integrate_chebychev_try(Expr* f, Expr* x);

/* `Integrate`ChebychevAlgebraic[f, x]` builtin.  Same algorithm as the cascade
 * entry.  Strict: returns NULL (unevaluated) when the integrand is not a
 * Chebychev binomial or the integral is non-elementary. */
Expr* builtin_integrate_chebychev(Expr* res);

/* Register the package builtin + attributes + docstring.  Called from
 * integrate_init(). */
void integrate_chebychev_init(void);

#endif /* INTEGRATE_CHEBYCHEV_H */
