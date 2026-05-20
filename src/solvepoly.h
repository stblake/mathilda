/*
 * solvepoly.h
 *
 * Single-variable polynomial-equality solver.  Backs the public
 * `Solve` router (src/solve.c) and is itself exposed as the
 * context-qualified builtin `Solve`SolvePolynomialEquality`.
 *
 * The solver follows the algorithm laid out in plans/SOLVE_PLAN.md:
 *   1. Move the equation to one side: poly = lhs − rhs.
 *   2. Recognise fast-path shapes (linear, quadratic, binomial,
 *      n-quadratic) and return closed-form rules.
 *   3. Otherwise: Collect → FactorSquareFree → Factor per square-free
 *      factor.  Each irreducible factor is handled per-degree:
 *        deg 1 → linear rule
 *        deg 2 → quadratic rules (discriminant-aware in Reals)
 *        binomial → all n radical roots
 *        n-quadratic → 2n radical roots via u = x^n
 *        deg 3 → Cardano (if Cubics → True) or 3 Root[] objects
 *        deg 4 → Ferrari (if Quartics → True) or 4 Root[] objects
 *        deg ≥ 5 → deg Root[] objects
 *      Multiplicities (from FactorSquareFree × Factor exponents) are
 *      preserved: roots are emitted once per unit of multiplicity.
 */

#ifndef SOLVEPOLY_H
#define SOLVEPOLY_H

#include <stdbool.h>
#include "expr.h"

typedef struct {
    bool cubics_radical;    /* default: false → emit Root[] for cubics  */
    bool quartics_radical;  /* default: false → emit Root[] for quartics */
} SolvePolyOpts;

/* Internal C-callable entry point.  Borrowed args: the caller retains
 * ownership of `equation`, `var`, and `dom`.  Returns a freshly-owned
 * solution `List[List[Rule[var, val]] ...]` on success, or NULL when
 * the input is not a polynomial equation in `var` (in which case the
 * caller should leave its enclosing Solve unevaluated). */
Expr* solvepoly_solve_polynomial_equality(Expr* equation,
                                          Expr* var,
                                          Expr* dom,
                                          const SolvePolyOpts* opts);

/* Built-in entry point for `Solve`SolvePolynomialEquality[...]`. */
Expr* builtin_solve_polynomial_equality(Expr* res);

/* Register the qualified builtin in the symbol table.  Idempotent. */
void solvepoly_init(void);

#endif /* SOLVEPOLY_H */
