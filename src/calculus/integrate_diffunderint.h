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
