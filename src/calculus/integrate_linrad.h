/* integrate_linrad.h
 *
 * Integration of rational functions of x and radicals of a common *linear*
 * argument -- the "linear radicals" (Euler-style) substitution.  Recognises an
 * integrand of the form
 *
 *     R[x, (a x + b)^(m1/n1), (a x + b)^(m2/n2), ...]
 *
 * where R is rational in its arguments, a, b are free of x, and every radical
 * shares the same linear base (a x + b).  With n = LCM[n1, n2, ...] the
 * substitution u = (a x + b)^(1/n) (so x = (u^n - b)/a, dx = (n/a) u^(n-1) du)
 * rationalises the integrand:
 *
 *     Integrate[R[...], x]
 *       = (n/a) Integrate[ R[(u^n - b)/a, u^M1, u^M2, ...] u^(n-1), u ],
 *         Mk = mk n / nk,
 *
 * a *rational function of u* that the BronsteinRational stage integrates
 * exactly.  The antiderivative is then back-substituted u -> (a x + b)^(1/n)
 * and verified by differentiation.  Examples the other stages miss cleanly:
 *
 *     Integrate[1/Sqrt[x + 1], x]      -> 2 Sqrt[1 + x]
 *     Integrate[1/(1 + x^(1/3)), x]    -> rational-in-u antiderivative
 *     Integrate[Sqrt[x]/(1 + Sqrt[x]), x]
 */
#ifndef INTEGRATE_LINRAD_H
#define INTEGRATE_LINRAD_H

#include "expr.h"

/* Direct entry point for the Integrate dispatcher's Automatic cascade.  Makes
 * the rationalising substitution and recurses into the full integrator.
 * Returns a freshly owned antiderivative of `f` with respect to symbol `x`,
 * or NULL when `f` is not a rational function of a single linear radical (or
 * the reduced rational integral does not close).  Does NOT take ownership of
 * `f` or `x`. */
Expr* integrate_linrad_try(Expr* f, Expr* x);

/* `Integrate`LinearRadicals[f, x]` builtin.  Same algorithm as the cascade
 * entry.  Strict: returns NULL (unevaluated) when the substitution does not
 * apply or the result fails to differentiate back to f. */
Expr* builtin_integrate_linrad(Expr* res);

/* Register the package builtin + attributes + docstring.  Called from
 * integrate_init(). */
void integrate_linrad_init(void);

#endif /* INTEGRATE_LINRAD_H */
