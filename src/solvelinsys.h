/*
 * solvelinsys.h
 *
 * Linear-system specialist for `Solve`.  Solves m linear equations in
 * n variables over Integers / Reals / Complexes via Gauss--Jordan
 * elimination on a reversed-column augmented matrix.  Backs `Solve`
 * for multi-variable inputs (And / List of Equals, or single Equal
 * with multiple unknowns) and is itself exposed as the context-
 * qualified builtin `Solve`SolveLinearSystem`.
 *
 * Output shape, matching Mathematica:
 *   - Unique solution           { { var1 -> val1, var2 -> val2, ... } }
 *   - Inconsistent system       { }
 *   - Tautology (no equations)  { { } }
 *   - Underdetermined system    { { pivot_vars -> expr in free vars } }
 *     (emits `Solve::svars` warning -- free vars produce no rule)
 *
 * Reversed-column ordering is what gives Mathematica's convention for
 * the under-determined case: with `vars = {x, y}` and one equation,
 * the rightmost listed variable (y) becomes the pivot and is solved
 * in terms of the leftmost (x).  See plans/SOLVE_PLAN.md for the algorithm
 * write-up.
 *
 * Memory contract is the same as solvepoly.c: every helper returns a
 * freshly-owned Expr* (or Expr**); inputs are borrowed and deep-copied
 * into the output.  NULL means "this specialist is inapplicable" and
 * the caller leaves its enclosing Solve unevaluated.
 */

#ifndef SOLVELINSYS_H
#define SOLVELINSYS_H

#include "expr.h"

/* Internal C-callable entry point.  Borrowed args -- the caller retains
 * ownership of `equations`, `vars`, and `dom`.
 *
 *   equations:  Equal[lhs, rhs], And[Equal[...], ...], or List[Equal[...], ...]
 *   vars:       List of symbols
 *   dom:        NULL (= Complexes), or one of Integers / Reals / Complexes
 *
 * Returns a freshly-owned solution list on success, or NULL when the
 * specialist is inapplicable (non-linear, malformed shape, or
 * unsupported domain). */
Expr* solvelinsys_solve_linear_system(Expr* equations,
                                      Expr* vars,
                                      Expr* dom);

/* Builtin entry: `Solve`SolveLinearSystem[eqns, vars]` and
 * `Solve`SolveLinearSystem[eqns, vars, dom]`. */
Expr* builtin_solve_linear_system(Expr* res);

/* Register the qualified builtin in the symbol table.  Idempotent. */
void solvelinsys_init(void);

#endif /* SOLVELINSYS_H */
