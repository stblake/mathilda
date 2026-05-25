/*
 * solverad.h
 *
 * Single-variable radical-equation solver -- the specialist dispatched
 * from `Solve` (src/solve.c) when the equation contains fractional-power
 * subterms (Sqrt, x^(p/q), nested radicals).  Also reachable directly as
 * the context-qualified builtin `Solve`SolveRadicalsEquality`.
 *
 * Algorithm:
 *
 *   1.  Compute e = Numerator(Together(lhs - rhs)).
 *   2.  Repeatedly locate radical atoms Power[base, p/q] (q > 1) anywhere
 *       in the working system (main equation + accumulated side
 *       equations).  For each distinct base g_i, introduce a fresh
 *       generator u_i = g_i^(1/L_i), where L_i = lcm of denominators of
 *       all exponents of g_i.  Replace every g_i^(p/q) with u_i^(p*L_i/q)
 *       in every equation, and append the side equation u_i^L_i - g_i = 0.
 *       Nested radicals are handled automatically by the loop -- a fresh
 *       atom inside a previously substituted base gets its own u_j in the
 *       next iteration.
 *   3.  Once no radicals remain anywhere, eliminate u_1, u_2, ... from
 *       the main equation by chained Resultant_{u_i}(main, side_eq_i, u_i).
 *       The chained order is the introduction order, so each side
 *       equation contributes exactly one fresh generator.
 *   4.  The eliminated main equation is now a polynomial in `var`; hand
 *       it off to Solve`SolvePolynomialEquality.
 *   5.  Verify every candidate by back-substitution into the *original*
 *       equation.  Symbolic 0 (Simplify) or near-zero N[] keeps the
 *       candidate; clearly non-zero N[] drops it; inconclusive results
 *       (free parameters) keep the candidate and raise Solve::nongen.
 *
 * Memory contract: caller retains ownership of `equation`, `var`, `dom`.
 * Returns a freshly-owned solution `List[List[Rule[var, val]] ...]` on
 * success, or NULL when the input has no radical atom in `var` (the
 * polynomial specialist is the correct dispatch target in that case).
 */

#ifndef SOLVERAD_H
#define SOLVERAD_H

#include "expr.h"

/* Internal C-callable entry point.  Borrowed args -- the caller retains
 * ownership of `equation`, `var`, and `dom`.
 *
 *   equation:  Equal[lhs, rhs]
 *   var:       a single Symbol
 *   dom:       NULL (= Complexes), or one of Reals / Complexes / Integers
 *
 * Returns a freshly-owned solution list on success, or NULL when the
 * specialist is inapplicable (no radicals, malformed shape, or
 * resultant elimination failed). */
Expr* solverad_solve_radicals_equality(Expr* equation,
                                       Expr* var,
                                       Expr* dom);

/* Builtin entry: `Solve`SolveRadicalsEquality[eq, var]` and
 * `Solve`SolveRadicalsEquality[eq, var, dom]`. */
Expr* builtin_solve_radicals_equality(Expr* res);

/* Register the qualified builtin in the symbol table.  Idempotent. */
void solverad_init(void);

#endif /* SOLVERAD_H */
