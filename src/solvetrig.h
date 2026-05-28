/*
 * solvetrig.h
 *
 * Trig canonicalisation pre-pass for the public `Solve` router
 * (src/solve.c).  Handles single-variable equations of the form
 * `f(Sin[var], Cos[var], Tan[var], ...) == 0` (and the hyperbolic
 * analogues) that the inverse-function specialist cannot peel
 * directly because more than one trig head over `var` appears.
 *
 * Algorithm:
 *
 *   1. residual = lhs - rhs.
 *   2. Apply TrigToExp to convert every trig head over `var` to a
 *      sum/product of Power[E, m * I * var] terms.
 *   3. Run Together to clear denominators, then take the numerator.
 *   4. Walk the numerator: replace every Power[E, m * I * var] with
 *      Power[var, m] (reusing `var` as the temporary u-symbol).
 *      Abort if any `var` occurrence is *not* inside such a Power,
 *      or if `m` is non-integer.
 *   5. Multiply through by `var ^ (-min m)` to clear any negative
 *      powers, then Expand.
 *   6. Solve as an ordinary polynomial equality in `var` (treated as
 *      `u`).  Roots are the candidate u-values.
 *   7. For each u-value `u_i`, build `Equal[Power[E, I*var], u_i]`
 *      and dispatch through the inverse-function specialist; that
 *      path produces the periodic family
 *          var = -I Log[u_i] + 2 Pi C[k].
 *   8. Aggregate.
 *
 * Mirrors Maxima's `trig-cannon` properties + `trig-subst-p`
 * pre-pass (Maxima rewrites to sin/cos; we route through Exp because
 * Mathilda already ships TrigToExp and the inverse-function
 * specialist's peel_exp).
 *
 * Memory contract: the caller retains ownership of every argument.
 * Returns a freshly-owned `List[List[Rule[var, val]], ...]` on
 * success, or NULL when the equation does not fit the trig-exp
 * shape (in which case the caller falls through to solverad).
 */

#ifndef SOLVETRIG_H
#define SOLVETRIG_H

#include <stdbool.h>
#include "expr.h"
#include "solveinv.h"  /* for SolveInvOpts */

/* Internal C-callable entry.  Borrowed args. */
Expr* solvetrig_solve_trig_equality(Expr* equation,
                                    Expr* var,
                                    Expr* dom,
                                    const SolveInvOpts* opts);

/* Cheap structural probe: returns true iff `expr` contains at least
 * one trig / hyperbolic head over `var` somewhere in its tree.
 * Used as a fast-fail guard so we never pay the TrigToExp /
 * Together / Expand cost on polynomial-only inputs. */
bool solvetrig_has_trig(const Expr* expr, const Expr* var);

/* Registration hook. */
void solvetrig_init(void);

#endif /* SOLVETRIG_H */
