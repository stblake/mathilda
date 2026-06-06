/* integrate_quadrad.h
 *
 * Integration of rational functions of x and radicals of a common *quadratic*
 * argument -- the classical Euler substitution.  Recognises an integrand of the
 * form
 *
 *     R[x, (a x^2 + b x + c)^(n1/2), (a x^2 + b x + c)^(n2/2), ...]
 *
 * where R is rational in its arguments, a, b, c are free of x, every radical is
 * a square root (denominator 2) and they share the same degree-2 radicand.  An
 * Euler substitution rationalises the integrand into a rational function of a
 * fresh variable u, which the BronsteinRational stage integrates exactly; the
 * antiderivative is then back-substituted.
 *
 * To keep results real-valued, exactly ONE substitution is chosen by the sign
 * of the leading coefficient a (and, when a < 0, the discriminant b^2 - 4 a c):
 *
 *     a > 0                         Euler I:   sqrt(a x^2+b x+c) = sqrt(a) x + u
 *     a < 0 and b^2 - 4 a c > 0     Euler III: sqrt(a x^2+b x+c) = (x - alpha) u
 *     a symbolic                    Euler I    (best-effort, a > 0 branch)
 *
 * Examples the other stages miss cleanly:
 *
 *     Integrate[1/Sqrt[x^2 + 1], x]    -> Log[x + Sqrt[1 + x^2]]
 *     Integrate[1/Sqrt[1 - x^2], x]    -> a real arcsine-equivalent
 *     Integrate[Sqrt[x^2 + x + 1], x]
 */
#ifndef INTEGRATE_QUADRAD_H
#define INTEGRATE_QUADRAD_H

#include "expr.h"

/* Direct entry point for the Integrate dispatcher's Automatic cascade.  Chooses
 * and applies a real-valued Euler substitution and recurses into the full
 * integrator.  Returns a freshly owned antiderivative of `f` with respect to
 * symbol `x`, or NULL when `f` is not a rational function of square roots of a
 * single quadratic argument (or the reduced rational integral does not close).
 * Does NOT take ownership of `f` or `x`. */
Expr* integrate_quadrad_try(Expr* f, Expr* x);

/* `Integrate`QuadraticRadicals[f, x]` builtin.  Same algorithm as the cascade
 * entry.  Strict: returns NULL (unevaluated) when the substitution does not
 * apply or the reduced rational integral does not close. */
Expr* builtin_integrate_quadrad(Expr* res);

/* Register the package builtin + attributes + docstring.  Called from
 * integrate_init(). */
void integrate_quadrad_init(void);

#endif /* INTEGRATE_QUADRAD_H */
