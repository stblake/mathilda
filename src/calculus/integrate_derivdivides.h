/* integrate_derivdivides.h
 *
 * Integration by substitution -- the classical "derivative-divides" method
 * (Maxima's `diffdiv`, Moses' SIN program).  Recognises an integrand of the
 * shape
 *
 *     f(x) = c . h(u(x)) . u'(x)
 *
 * for some kernel u(x), reduces the problem to  Integrate[h[u], u]  and
 * back-substitutes  u -> u(x).  Examples that the rational / Risch / table
 * stages miss:
 *
 *     Integrate[Sin[x] Sqrt[1 - Cos[x]], x]   -> 2/3 (1 - Cos[x])^(3/2)
 *     Integrate[x Exp[x^2], x]                -> Exp[x^2]/2
 *     Integrate[1/(x Log[x]), x]              -> Log[Log[x]]
 *
 * Two complementary strategies, each followed by an unconditional
 * verification-by-differentiation gate (Simplify[D[result, x] - f] === 0),
 * which is what selects the correct branch when a substitution introduces
 * multiple inverse-function branches:
 *
 *   (1) Direct quotient:  q = Cancel[Together[f / D[u(x), x]]], then
 *       q /. u(x) -> u; accept when the result is free of x.  Picks the
 *       correct radical branch inherently (no squaring) and handles
 *       transcendental compositions such as Exp[x^2].  Quiet and cheap.
 *
 *   (2) Eliminate / Solve:  build the differential relation via
 *       Eliminate[{Dt[y] == f Dt[x], u == u(x), Dt[u == u(x)]}, {x, Dt[x]}],
 *       Solve for Dt[y], reduce each branch with Factor //@ , PowerExpand and
 *       Cancel[. / Dt[u]], and integrate the branch that differentiates back
 *       to f.  Handles algebraic entanglement the direct method cannot (e.g.
 *       1/Cos[x] canonicalised to Sec[x], or the Pythagorean reduction of
 *       Sin[x]^3 over u = Cos[x]).  May emit Eliminate's ::ifun / ::alg
 *       branch-caveat diagnostics, so it is reserved for the explicit method.
 */
#ifndef INTEGRATE_DERIVDIVIDES_H
#define INTEGRATE_DERIVDIVIDES_H

#include "expr.h"

/* Direct entry point for the Integrate dispatcher's Automatic cascade.
 * Runs the *direct quotient* strategy only (strategy 1 above): fast, emits
 * no diagnostics, and picks the correct branch inherently.  Returns a
 * freshly owned antiderivative of `f` with respect to symbol `x`, or NULL
 * when no kernel substitution closes the integral.  Does NOT take ownership
 * of `f` or `x`. */
Expr* integrate_derivdivides_try(Expr* f, Expr* x);

/* `Integrate`DerivativeDivides[f, x]` builtin.  Runs the direct strategy and
 * then the more thorough Eliminate/Solve branch-search (strategy 2).  Strict:
 * returns NULL (unevaluated) when neither strategy produces a result that
 * differentiates back to f. */
Expr* builtin_integrate_derivdivides(Expr* res);

/* Register the package builtin + attributes + docstring.  Called from
 * integrate_init(). */
void integrate_derivdivides_init(void);

#endif /* INTEGRATE_DERIVDIVIDES_H */
