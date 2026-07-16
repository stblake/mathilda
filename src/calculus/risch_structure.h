/* risch_structure.h — the Risch structure theorems (Bronstein §9.3).
 *
 * The structure theorems (Bronstein, *Symbolic Integration I*, 2nd ed., §9.3,
 * Corollaries 9.3.1 (complex) / 9.3.2 (real)) reduce the decision problems that
 * gate every tower extension to a Q-LINEAR SYSTEM over the constants:
 *
 *   - Is a proposed new logarithm log(a) algebraically independent of the tower,
 *     or reducible?  Reducible iff Da/a is a Q-linear combination of the
 *     monomials' generators (Cor. 9.3.1(i), eq. 9.8).
 *   - Is a proposed new exponential exp(b) independent, or reducible?  Reducible
 *     iff Db is a Q-linear combination of the generators (Cor. 9.3.1(ii),
 *     eq. 9.9) — equivalently, Db is the logarithmic derivative of a K-radical.
 *
 * The "generator" of a monomial t_i is Dt_i for a logarithm and Dt_i/t_i for an
 * exponential (both elements of the field below t_i); the target is Da/a
 * (log test) or Db (exp test).  This module implements the underlying rational
 * Q-span decision and the two structure-theorem front-ends.
 *
 * Builtins (Q-span decision core + structure-theorem tests):
 *   Risch`RationalSpan[theta, {g1,...,gm}, {vars}]
 *        -> {r1,...,rm} (rational coefficients with theta = Sum r_i g_i) or False.
 *   Risch`LogReducible[a, x, {{t_i, "Exp"|"Log", Dt_i}, ...}]
 *        -> the Q-coefficients if log(a) is reducible over the tower, else False.
 *   Risch`ExpReducible[b, x, {{t_i, "Exp"|"Log", Dt_i}, ...}]
 *        -> the Q-coefficients if exp(b) is reducible over the tower, else False.
 *
 * Memory contract: the standard Mathilda BuiltinFunc rule.
 */

#ifndef MATHILDA_RISCH_STRUCTURE_H
#define MATHILDA_RISCH_STRUCTURE_H

#include <stdbool.h>
#include <stddef.h>
#include "expr.h"

/* Decide whether theta = Sum_i r_i gens[i] for rational r_i, matching
 * coefficients over the variables `vars` (a List Expr).  Returns an owned
 * List[r_0,...,r_{m-1}] of rational numbers on success, or NULL if theta is not
 * in the rational span of the generators. */
Expr* risch_rational_span(const Expr* theta, Expr* const* gens, size_t m,
                          const Expr* vars);

/* Register the Risch` structure-theorem builtins.  Called from integrate_init(). */
void integrate_risch_structure_init(void);

#endif /* MATHILDA_RISCH_STRUCTURE_H */
