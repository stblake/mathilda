/*
 * integrate_newton_leibniz.h -- Definite integration by the Newton-Leibniz
 * rule (the fundamental theorem of calculus).
 *
 * Given an integrand f, a variable x, and bounds a, b, this module:
 *   1. computes an antiderivative F = Integrate[f, x] with the existing
 *      indefinite cascade,
 *   2. locates the real singular points of f and F strictly inside (a, b)
 *      (poles of the rational part -- see the detector below),
 *   3. splits [a, b] at those points and forms the telescoping sum
 *          Sum_i ( F(p_{i+1}^-) - F(p_i^+) ),
 *      evaluating each boundary with the Limit engine (one-sided at
 *      singular / infinite points) so improper integrals get their correct
 *      value and genuinely divergent integrals report Infinity /
 *      ComplexInfinity / Indeterminate.
 *
 * The definite form is reachable three ways:
 *   Integrate[f, {x, a, b}]                       (automatic dispatch)
 *   Integrate[f, {x, a, b}, Method -> "NewtonLeibniz"]
 *   Integrate`NewtonLeibniz[f, {x, a, b}]         (explicit entry point)
 *
 * The pole detector is exposed for inspection / unit testing / reuse:
 *   Integrate`SingularPoints[expr, {x, a, b}]     -> sorted List of the
 *                                                    interior real poles.
 */

#ifndef MATHILDA_INTEGRATE_NEWTON_LEIBNIZ_H
#define MATHILDA_INTEGRATE_NEWTON_LEIBNIZ_H

#include "expr.h"

/*
 * Core entry point for the Integrate dispatcher's definite path.
 *   f, x, a, b are borrowed (not consumed).
 *   `method` is a NUL-terminated indefinite-method name to pass through to
 *   the inner Integrate[f, x, Method -> method] call, or NULL for the default
 *   Automatic cascade.
 * Returns a freshly-allocated Expr* (the definite value, possibly a
 * divergence symbol), or NULL to leave the definite integral unevaluated
 * (unknown antiderivative, undecidable pole position, or a boundary limit
 * the engine could not determine).
 */
Expr* integrate_newton_leibniz_try(Expr* f, Expr* x, Expr* a, Expr* b,
                                   const char* method);

/* Emit the Mathematica-style `Integrate::idiv: Integral of <f> does not
 * converge on {<a>, <b>}.` warning to stderr.  Shared so other definite
 * mechanisms (e.g. the residue method) report a conclusively divergent
 * integral with the same diagnostic. */
void integrate_emit_idiv(Expr* f, Expr* a, Expr* b);

/* `Integrate`NewtonLeibniz[f, {x, a, b}]` builtin.  Strict: returns NULL on
 * any non-applicable input. */
Expr* builtin_integrate_newton_leibniz(Expr* res);

/* `Integrate`SingularPoints[expr, {x, a, b}]` builtin.  Returns the sorted
 * List of confirmed real poles of `expr` strictly inside (a, b), or NULL on
 * malformed input. */
Expr* builtin_integrate_singular_points(Expr* res);

/* Register the package builtins + attributes + docstrings. */
void integrate_newton_leibniz_init(void);

#endif /* MATHILDA_INTEGRATE_NEWTON_LEIBNIZ_H */
