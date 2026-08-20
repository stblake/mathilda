/*
 * solvealways.h
 *
 * `SolveAlways[eqns, vars]` -- find values of parameters (the symbols
 * appearing in `eqns` but not in `vars`) that make every equation hold
 * for all values of `vars`.
 *
 * Reduction: each `lhs == rhs` is rewritten as the polynomial `lhs - rhs`
 * in `vars`; the coefficients of that polynomial must all vanish for the
 * equation to be an identity in `vars`.  Collect every such coefficient,
 * then hand the system to `Solve` with the parameters as unknowns.
 */

#ifndef SOLVEALWAYS_H
#define SOLVEALWAYS_H

#include "expr.h"

Expr* builtin_solvealways(Expr* res);

/* Register `SolveAlways` and its Protected attribute.  Idempotent. */
void solvealways_init(void);

#endif /* SOLVEALWAYS_H */
