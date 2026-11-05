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
#include <stdbool.h>

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

/* As above, but when `principal_value` is true and the integrand has interior
 * real poles, computes the Cauchy principal value: valid only when every
 * interior pole is of odd order (the integrand changes sign across it, so the
 * one-sided divergences cancel), in which case PV = Re[F(b) - F(a)] with the
 * real-branch antiderivative.  An even-order interior pole has no principal
 * value and yields NULL (Integrate::idiv).  With no interior pole this is just
 * the ordinary definite integral.
 *
 * `assumptions` (may be NULL) is threaded into pole classification (so a
 * symbolic root like x = a against a symbolic bound is decided by Refine rather
 * than failing undecidable) and into the boundary evaluation (Limit and a Refine
 * of the substituted endpoint), which is what lets a symbolic-endpoint integral
 * such as Integrate[1/Sqrt[a^2-x^2], {x,0,a}, Assumptions -> a>0] close. */
Expr* integrate_newton_leibniz_try_pv(Expr* f, Expr* x, Expr* a, Expr* b,
                                      const char* method, bool principal_value,
                                      Expr* assumptions);

/* True iff `e`, with its free variables instantiated at a sample value, N[.]'s
 * to a finite complex number with a nonzero imaginary part.  A non-numericizing
 * (genuinely symbolic) or real result returns false.  Used by the definite
 * dispatcher so a Refine cleanup does not turn a real result (e.g.
 * Re[CosIntegral[a]-CosIntegral[-a]]) into a spurious complex one. */
bool nl_numeric_is_complex(Expr* e);

/* Emit the Mathematica-style `Integrate::idiv: Integral of <f> does not
 * converge on {<a>, <b>}.` warning to stderr.  Shared so other definite
 * mechanisms (e.g. the residue method) report a conclusively divergent
 * integral with the same diagnostic. */
void integrate_emit_idiv(Expr* f, Expr* a, Expr* b);

/* Hold the idiv diagnostic of one cascade region.  Under PrincipalValue the
 * Newton-Leibniz stage reports idiv when it cannot form the PV, but a later
 * stage (Roach's PV theorem) may still close the integral, so the definite
 * dispatcher brackets that stage with begin/end: an idiv emitted inside the
 * region (by the stage or a sub-integral) is held, not printed, and end()
 * says whether one was.  The dispatcher re-emits it only if nothing answered.
 * Regions nest (begin returns the enclosing region's flag, end restores it);
 * the nesting counter is registered with TimeConstrained, so a timeout that
 * unwinds through a region does not leave idiv held for the session. */
bool integrate_idiv_defer_begin(void);
bool integrate_idiv_defer_end(bool prev);
int* integrate_idiv_defer_counter(void);

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
