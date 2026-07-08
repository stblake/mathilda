/*
 * integrate_ramanujan.h -- Definite integration of half-line integrals
 * Integrate[f, {x, 0, Infinity}] by the Mellin-transform / Ramanujan Master
 * Theorem method.
 *
 * Ramanujan's Master Theorem: for a function with the Maclaurin expansion
 *   f(x) = Sum_{k>=0} (-1)^k phi(k) x^k / k!,
 * the Mellin transform is
 *   Integrate[x^(s-1) f(x), {x, 0, Infinity}] = Gamma(s) phi(-s),
 * valid on the fundamental strip 0 < Re(s) < strip under Hardy's conditions.
 *
 * This module realizes the theorem in two layers:
 *
 *   Layer 1 (table + operational rules): a small table of proven base Mellin
 *     transforms M_f(s) = Integrate[x^(s-1) f(x), {x,0,Infinity}] (exponential,
 *     Gaussian, algebraic binomial (p + q x^m)^(-a), Cos, Sin, BesselJ), each
 *     carrying its convergence strip.  A definite integrand is decomposed
 *     term-by-term into  C * x^rho * kernel(x); the kernel is matched against
 *     the table, the power prefactor shifts s = rho + 1, and internal scalings
 *     fold into the transform.  Every application is gated on its strip, so the
 *     result is unconditionally correct.
 *
 *   Layer 2 (guarded general RMT): reserved for integrands whose Maclaurin
 *     coefficient function is a hypergeometric-type ratio (the class where
 *     Hardy's conditions provably hold).  Currently conservative -- Layer 1
 *     handles the shipped families; Layer 2 is a documented extension point.
 *
 * Verification is SYMBOLIC and correct-by-construction: each table identity is
 * a theorem, and the strip gate (checked numerically for numeric s, or via the
 * supplied Assumptions for symbolic s) guarantees convergence.  There is NO
 * numeric (NIntegrate) crosscheck anywhere in this path (project rule).  When
 * the strip cannot be established, or the integrand is out of scope (e.g. a
 * product of two transcendental kernels -- Mellin convolution / Meijer-G
 * territory), the integral is left unevaluated (NULL).
 *
 * Reachable three ways:
 *   Integrate[f, {x, 0, Infinity}]                                  (automatic
 *                                       cascade, after residue + Newton-Leibniz)
 *   Integrate[f, {x, 0, Infinity}, Method -> "RamanujanMasterTheorem"]
 *   Integrate`RamanujanMasterTheorem[f, {x, 0, Infinity}]           (explicit)
 */

#ifndef MATHILDA_INTEGRATE_RAMANUJAN_H
#define MATHILDA_INTEGRATE_RAMANUJAN_H

#include "expr.h"

/*
 * Core entry point for the Integrate dispatcher's definite path.
 *   f, x, a, b, assumptions are borrowed (not consumed); assumptions may be
 *   NULL.  x must be a symbol.  Applies only when a is the literal 0 and b is
 *   +Infinity; returns NULL for every other bound.
 * Returns a freshly-allocated Expr* (the closed-form definite value), or NULL
 * to leave the definite integral unevaluated.
 */
Expr* integrate_ramanujan_try(Expr* f, Expr* x, Expr* a, Expr* b,
                              Expr* assumptions);

/* `Integrate`RamanujanMasterTheorem[f, {x,0,Infinity}]` (optionally with
 * Assumptions -> ...) builtin.  Strict: returns NULL on any non-applicable
 * input. */
Expr* builtin_integrate_ramanujan(Expr* res);

/* Register the package builtin + attributes + docstring. */
void integrate_ramanujan_init(void);

#endif /* MATHILDA_INTEGRATE_RAMANUJAN_H */
