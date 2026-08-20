/*
 * solveint.h
 *
 * Integer-domain (Diophantine) solving pre-pass for
 * Solve[eqns && constraints, vars, Integers].
 *
 * Contract (mirrors solvemod.h):
 *   solveint_solve_integer(expr, vars, dom) returns a canonical solution
 *   List[List[Rule...], ...] on success (ascending by value tuple),
 *   the empty List {} for a provably empty solution set, or NULL to
 *   decline -- in which case Solve stays unevaluated / falls through to
 *   the ordinary dispatch.  All arguments are borrowed; the returned
 *   Expr* is freshly owned by the caller.
 *
 *   `dom` must be the symbol Integers for this pre-pass to engage; the
 *   caller (builtin_solve) guarantees that before calling.
 *
 * Phase 1 scope: equations that are polynomial (after clearing rational
 * denominators) in the vars, together with inequality / ordering /
 * disequation constraints that bound every variable to a finite integer
 * range.  Solved by bound propagation + recursive elimination with an
 * exact univariate leaf.  Unbounded or non-polynomial inputs are
 * declined (NULL).  Linear-Diophantine parametric output, Pell, and the
 * research-grade forms are later phases.
 */

#ifndef SOLVEINT_H
#define SOLVEINT_H

#include "expr.h"

/* Top-level integer-domain solver.  See the contract above. */
Expr* solveint_solve_integer(Expr* expr, Expr* vars, Expr* dom);

/* Registers Solve`SolveIntegers (the independently-testable entry). */
void solveint_init(void);

#endif /* SOLVEINT_H */
