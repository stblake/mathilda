/*
 * solve.h
 *
 * `Solve` -- the public equation-solving builtin.  Acts as a router
 * that parses options, normalises arguments, and dispatches to a
 * specialist (currently Solve`SolvePolynomialEquality in solvepoly.c).
 *
 * Architecture: docs/notes/plans/SOLVE_PLAN.md.
 */

#ifndef SOLVE_H
#define SOLVE_H

#include "expr.h"

Expr* builtin_solve(Expr* res);

/* Register Solve, the option-name docstrings, and the polynomial-
 * equality specialist.  Idempotent. */
void solve_init(void);

/* Solve#6: true when this Solve call may report Solve::ifun (once per
 * top-level Solve); marks it reported. */
bool solve_ifun_claim(void);

#endif /* SOLVE_H */
