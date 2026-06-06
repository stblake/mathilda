/* integrate_linratiorad.h
 *
 * Integration of rational functions of x and radicals of a common *linear
 * fractional* (Möbius) argument -- the "linear ratio radicals" rationalising
 * substitution.  Recognises an integrand of the form
 *
 *     R[x, ((a x + b)/(c x + d))^(n1/m1), ((a x + b)/(c x + d))^(n2/m2), ...]
 *
 * where R is rational in its arguments, a, b, c, d are free of x with
 * a d - b c != 0, and every radical shares the same linear-fractional base
 * r = (a x + b)/(c x + d).  With m = LCM[m1, m2, ...] the substitution
 * u = r^(1/m) (so r = u^m) inverts the Möbius map to
 *
 *     x  = (d u^m - b)/(a - c u^m),
 *     dx = m (a d - b c) u^(m-1)/(a - c u^m)^2 du,
 *
 * which rationalises the integrand:
 *
 *     Integrate[R[...], x]
 *       = Integrate[ R[(d u^m - b)/(a - c u^m), u^M1, ...]
 *                    · m (a d - b c) u^(m-1)/(a - c u^m)^2, u ],  Mk = nk m/mk,
 *
 * a *rational function of u* that the BronsteinRational stage integrates
 * exactly.  The antiderivative is then back-substituted
 * u -> ((a x + b)/(c x + d))^(1/m).  The Möbius substitution is an exact
 * bijection on the relevant domain, so the result is correct by construction
 * once the rational sub-integral closes -- no differentiate-back check is
 * performed (matching integrate_linrad.c / integrate_quadrad.c).  Examples the
 * linear- and quadratic-radical stages miss cleanly:
 *
 *     Integrate[1/Sqrt[(x + 1)/(x - 1)], x]
 *     Integrate[((x - 1)/(x + 1))^(1/3)/x, x]
 */
#ifndef INTEGRATE_LINRATIORAD_H
#define INTEGRATE_LINRATIORAD_H

#include "expr.h"

/* Direct entry point for the Integrate dispatcher's Automatic cascade.  Makes
 * the rationalising substitution and recurses into the full integrator.
 * Returns a freshly owned antiderivative of `f` with respect to symbol `x`,
 * or NULL when `f` is not a rational function of a single linear-fractional
 * radical (or the reduced rational integral does not close).  Does NOT take
 * ownership of `f` or `x`. */
Expr* integrate_linratiorad_try(Expr* f, Expr* x);

/* `Integrate`LinearRatioRadicals[f, x]` builtin.  Same algorithm as the cascade
 * entry.  Strict: returns NULL (unevaluated) when the substitution does not
 * apply or the reduced rational integral does not close. */
Expr* builtin_integrate_linratiorad(Expr* res);

/* Register the package builtin + attributes + docstring.  Called from
 * integrate_init(). */
void integrate_linratiorad_init(void);

#endif /* INTEGRATE_LINRATIORAD_H */
