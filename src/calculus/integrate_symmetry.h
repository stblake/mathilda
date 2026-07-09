/*
 * integrate_symmetry.h -- Definite integration by symmetry reduction over an
 * interval symmetric about the origin, [-c, c] (c finite or +Infinity).
 *
 * Two reductions, each verified symbolically (Simplify-to-0, never
 * PossibleZeroQ) and each gated on the *half* integral actually converging so
 * a divergent principal value is never mistaken for a value:
 *
 *   odd  f(-x) == -f(x):  Integrate[f, {x, -c, c}] = 0
 *                         (only when Integrate[f, {x, 0, c}] is finite)
 *   even f(-x) ==  f(x):  Integrate[f, {x, -c, c}] = 2 Integrate[f, {x, 0, c}]
 *                         (finite intervals only; the infinite even half-line is
 *                          already handled by the residue even reduction)
 *
 * This is the first mechanism the definite dispatcher tries: an odd integrand
 * short-circuits to 0 without any antiderivative, and both reductions map the
 * problem onto the (usually easier) half-line / half-interval.  On NULL the
 * dispatcher falls through to residue / Newton-Leibniz / Ramanujan / DiffUnderInt
 * exactly as before, so it can only add coverage, never remove it.
 *
 * Reachable three ways:
 *   Integrate[f, {x, -c, c}]                        (automatic dispatch, first)
 *   Integrate[f, {x, -c, c}, Method -> "Symmetry"]  (strict: no fallback)
 *   Integrate`Symmetry[f, {x, -c, c}]               (explicit entry point)
 */

#ifndef MATHILDA_INTEGRATE_SYMMETRY_H
#define MATHILDA_INTEGRATE_SYMMETRY_H

#include "expr.h"

/*
 * Core entry for the definite dispatcher.  f, x, a, b, assumptions are borrowed
 * (not consumed; assumptions may be NULL).  Returns a freshly-allocated Expr*
 * (the definite value) on success, or NULL to leave the reduction to another
 * mechanism (interval not symmetric, integrand neither odd nor even, or the half
 * integral did not close to a finite value).
 */
Expr* integrate_symmetry_try(Expr* f, Expr* x, Expr* a, Expr* b,
                             Expr* assumptions);

/* `Integrate`Symmetry[f, {x, a, b}]` builtin.  Strict: NULL on any
 * non-applicable input. */
Expr* builtin_integrate_symmetry(Expr* res);

/* Register the builtin + attributes + docstring. */
void integrate_symmetry_init(void);

#endif /* MATHILDA_INTEGRATE_SYMMETRY_H */
