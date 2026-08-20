/*
 * solvenlsys.h
 *
 * Nonlinear polynomial-system specialist for `Solve` (src/solve.c).
 * Also reachable directly as the context-qualified builtin
 * `Solve`SolveNonlinearSystem`.
 *
 * Where solvelinsys.c handles affine systems by Gauss--Jordan, this
 * specialist handles genuinely NONLINEAR systems of polynomial
 * equations over Q for a ZERO-DIMENSIONAL solution set: it computes a
 * lexicographic Gröbner basis (which is triangular for a shape-position
 * ideal), solves the univariate generator in the last variable with the
 * single-variable polynomial solver (src/poly/solvepoly.c), and
 * back-substitutes upward, branching over each root and verifying every
 * completed tuple against the original equations.
 *
 * Output shape, matching Mathematica and solvelinsys.c:
 *   - Finite solutions          { {v1 -> a1, ...}, {v1 -> b1, ...}, ... }
 *   - Provably inconsistent      { }   (unit ideal, or no domain points)
 *   - Tautology (no equations)   { { } }
 *
 * Reserved / not handled (leave Solve unevaluated -> return NULL):
 *   - non-polynomial systems (a transcendental head or a symbol that is
 *     not one of `vars`),
 *   - positive-dimensional ideals (infinitely many solutions), which
 *     emit the advisory `Solve::nsdim`,
 *   - any branch the univariate solver cannot reduce (never emit a false
 *     `{}` -- an incomplete search returns NULL instead).
 *
 * Memory contract mirrors solvelinsys.c: `equations`, `vars`, `dom`, and
 * `opts` are BORROWED; every returned Expr* is freshly owned.  NULL means
 * "this specialist is inapplicable" and the caller leaves Solve
 * unevaluated.
 */

#ifndef SOLVENLSYS_H
#define SOLVENLSYS_H

#include "expr.h"
#include "solvepoly.h"   /* SolvePolyOpts (Cubics / Quartics passthrough) */

/* Internal C-callable entry point.  Borrowed args.
 *
 *   equations:  Equal[lhs, rhs], And[Equal[...], ...], or List[Equal[...], ...]
 *   vars:       List of >= 1 distinct symbols
 *   dom:        NULL (= Complexes), or one of Integers / Reals / Complexes
 *   opts:       Cubics / Quartics radical flags forwarded to the univariate
 *               step (may be NULL, in which case defaults are used)
 *
 * Returns a freshly-owned solution list on success, `{}`/`{{}}` for the
 * inconsistent/tautology cases, or NULL when the specialist is
 * inapplicable (non-polynomial, positive-dimensional, or an unsolved
 * branch). */
Expr* solvenlsys_solve_nonlinear_system(Expr* equations,
                                        Expr* vars,
                                        Expr* dom,
                                        const SolvePolyOpts* opts);

/* Builtin entry: `Solve`SolveNonlinearSystem[eqns, vars]` and
 * `Solve`SolveNonlinearSystem[eqns, vars, dom]`. */
Expr* builtin_solve_nonlinear_system(Expr* res);

/* Register the qualified builtin in the symbol table.  Idempotent. */
void solvenlsys_init(void);

#endif /* SOLVENLSYS_H */
