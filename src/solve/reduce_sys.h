/*
 * reduce_sys.h
 *
 * Parametric linear-system engine for `Reduce` (REDUCE_PLAN.md, Phase 4).
 *
 * Solves a system of LINEAR equations in the solve variables over the complex
 * numbers, where the coefficients may involve free parameters, and reports the
 * COMPLETE solution set with case analysis on the parameters -- the piece Solve
 * omits.  It runs symbolic Gaussian elimination: a nonzero-constant pivot is
 * used directly; a symbolic pivot `p` splits into `p != 0` (use it) and `p == 0`
 * (solve that parameter relation, substitute, and recurse).  Back-substitution
 * gives each variable as a function of the parameters.
 *
 *   Reduce[a x + y == 1 && x + y == 0, {x, y}]
 *       ->  a != 1 && x == 1/(a-1) && y == 1/(1-a)
 *   Reduce[a x == 1 && x == 2, {x}]
 *       ->  2 a - 1 == 0 && x == 2
 *
 * Declines (returns NULL) on a non-linear equation, an inequality/`!=` atom
 * (meaningless over Complexes), or a parameter relation it cannot solve
 * rationally.
 */
#ifndef REDUCE_SYS_H
#define REDUCE_SYS_H

#include "expr.h"
#include "reduce_form.h"

/* Solve the DNF `F` (each conjunct a linear equation system) in `vars` over the
 * complex numbers.  Returns a freshly-owned Expr (True / False / the complete
 * case-split solution, OR-ed over disjuncts), or NULL to decline. */
Expr* reduce_eq_system(const RForm* F, Expr** vars, int nv);

#endif /* REDUCE_SYS_H */
