/*
 * integrate_diffunderint.h -- Definite integration by differentiation under
 * the integral sign (the Leibniz rule / "Feynman's trick").
 *
 * For a definite integral I(p) = Integrate[f(x,p), {x,a,b}] that depends on a
 * free parameter p, this module:
 *   1. differentiates the integrand w.r.t. p (Leibniz rule),
 *   2. evaluates the (usually simpler) inner integral J(p) = Integrate[
 *      D[f,p], {x,a,b}] with the existing definite engine,
 *   3. integrates J(p) back over the parameter and fixes the constant of
 *      integration with an EXACT base value I(p0) (a value of p where the
 *      integral is known -- typically where f|_{p->p0} vanishes identically,
 *      or reduces to a directly-integrable form).
 *
 * Every case is a first-order ODE in the parameter, I'(p) = lambda(p) I(p) +
 * M(p) (Boulnois, arXiv:2308.09619).  The common degenerate case lambda = 0 is
 * pure quadrature (Stage A).  A genuine lambda != 0 (e.g. the Gaussian
 * Integrate[Exp[-x^2] Cos[2 a x], {x,0,Inf}], Conrad section 6) is recognized
 * via an integration-by-parts boundary identity and solved with an integrating
 * factor (Stage B).  Conditional (piecewise) inner integrals produce Piecewise
 * / Min / Max output (Stage C).  Specialized closers handle self-similar
 * scaling substitutions and the power-raising recurrence (Stage D).
 *
 * Three FINITE-DOMAIN families need neither a pre-existing parameter nor an
 * engine-safe inner integral, so they have dedicated closers that canonicalise
 * with a change of variables, INTRODUCE an artificial parameter, and supply the
 * inner integral in closed form (never routed through the general engine, which
 * hangs -- or, for the radical case, returns a WRONG value -- on them):
 *   Family 1 (power-log):      Log[1 + c x^p]/(x Sqrt[1-x^(2p)]) on {0,1}
 *                              -> (Pi^2/8 - ArcCos[c]^2/2)/p   [u = x^p]
 *   Family 2 (secant-radical): Sec[2x] Log[1 + c Sqrt[1-Tan[x]^2]] on {0,Pi/4}
 *                              -> Pi^2/8 - ArcCos[c]^2/2        [t = Tan[x]]
 *   Family 3 (tangent-power):  Csc[2x]^2 Log[1 + Tan[x]^a] on {0,Pi/4}
 *                              -> (Pi Csc[Pi/a] - a)/4          [t = Tan[x]]
 * Families 1 and 2 share the ArcCos[q]/Sqrt[1-q^2] inner integral; family 3 is a
 * direct Beta/digamma reflection with a rational anchor (a0 = 3) self-check.
 *
 * Verification is SYMBOLIC and correct-by-construction: PossibleZeroQ[
 * D[I,p] - J] plus an exact base value.  There is NO numeric (NIntegrate)
 * crosscheck anywhere in this code path (project rule).  The conditional-
 * convergence pitfall (Conrad section 12) is caught automatically: a
 * non-integrable D[f,p] makes the inner Integrate fail to close, so that
 * parameter is skipped.
 *
 * Reachable three ways:
 *   Integrate[f, {x,a,b}]                                 (automatic cascade,
 *                                                          after residue + FTC)
 *   Integrate[f, {x,a,b}, Method -> "DiffUnderInt"]
 *   Integrate`DiffUnderInt[f, {x,a,b}]                    (explicit entry point)
 */

#ifndef MATHILDA_INTEGRATE_DIFFUNDERINT_H
#define MATHILDA_INTEGRATE_DIFFUNDERINT_H

#include "expr.h"

/*
 * Core entry point for the Integrate dispatcher's definite path.
 *   f, x, a, b, assumptions are borrowed (not consumed); assumptions may be
 *   NULL.  x must be a symbol; a, b are the (possibly infinite) bounds.
 * Returns a freshly-allocated Expr* (the closed-form definite value), or NULL
 * to leave the definite integral unevaluated.
 */
Expr* integrate_diffunderint_try(Expr* f, Expr* x, Expr* a, Expr* b,
                                 Expr* assumptions);

/* `Integrate`DiffUnderInt[f, {x,a,b}]` (optionally with Assumptions -> ...)
 * builtin.  Strict: returns NULL on any non-applicable input. */
Expr* builtin_integrate_diffunderint(Expr* res);

/* Register the package builtin + attributes + docstring. */
void integrate_diffunderint_init(void);

#endif /* MATHILDA_INTEGRATE_DIFFUNDERINT_H */
