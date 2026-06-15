/* intsimp.h — simplification helpers for the rational-function
 * integrator (`Integrate`` package).
 *
 * Carved out of intrat.c so that the simplification surface — both
 * the internal Sqrt / sign / canonic-zero machinery used while
 * building the integral, and the post-pipeline cleanup applied to
 * the resulting integral — lives in one focused module.
 *
 * Two flavours of helper live here:
 *
 *   1. Algorithmic simplification used *inside* the integrator:
 *        intsimp_pos_sqrt          — simplify Sqrt[e] under the
 *                                    positive-symbol assumption
 *        intsimp_pos_sqrt_factor   — same, applied factor-wise
 *        intsimp_sign_pos_assumption / intsimp_numeric_sign
 *                                  — sign tests fed into the Sqrt
 *                                    simplifier and LogToReal
 *                                    dispatcher
 *        intsimp_zero_q            — Cancel[Together[e]] === 0
 *        intsimp_has_radical       — does e contain a Sqrt / Power[…,
 *                                    Rational[_,_]] subexpression?
 *        intsimp_simplify_if_radical
 *                                  — apply Simplify[] when the input
 *                                    carries an unsimplified radical
 *                                    (takes ownership)
 *
 *   2. Output cleanup applied to the resulting integral:
 *        intsimp_log_to_arctanh    — combine c Log[A] ± c Log[B]
 *                                    pairs into single Log / ArcTanh
 *        intsimp_strip_log_constants
 *                                  — Log[c·p] -> Log[p] when FreeQ[c,x]
 *        intsimp_distribute_plus   — expand Times-over-Plus
 *        intsimp_normalize_inverse_trig_signs
 *                                  — ArcTan[-arg] -> -ArcTan[arg]
 *                                    (and ArcTanh likewise)
 *
 * Memory contract: every helper returns a freshly-allocated Expr*
 * the caller owns; none of them free their input arguments unless
 * explicitly noted (intsimp_simplify_if_radical takes ownership).
 */

#ifndef MATHILDA_INTSIMP_H
#define MATHILDA_INTSIMP_H

#include <stdbool.h>
#include "expr.h"

/* --- Algorithmic simplification ----------------------------------- */

Expr* intsimp_pos_sqrt(Expr* e);
Expr* intsimp_pos_sqrt_factor(Expr* e);

int   intsimp_sign_pos_assumption(Expr* e);
int   intsimp_numeric_sign(Expr* e);

bool  intsimp_zero_q(Expr* e);

bool  intsimp_has_radical(const Expr* e);

/* Takes ownership of `e`.  If e contains a radical, returns
 * Simplify[e] (with e freed); otherwise returns e unchanged. */
Expr* intsimp_simplify_if_radical(Expr* e);

/* --- Output cleanup ----------------------------------------------- */

Expr* intsimp_log_to_arctanh(Expr* e, Expr* x);
Expr* intsimp_strip_log_constants(Expr* e, Expr* x);
Expr* intsimp_distribute_plus(Expr* e);
Expr* intsimp_normalize_inverse_trig_signs(Expr* e);

/* Final normalisation applied once to every closed antiderivative, scoped
 * to radical-bearing results (those carrying a Power[f(x), p/q] of the
 * integration variable).  For such results it flattens x-free
 * reciprocal/product powers ((1/a)^(1/3) -> a^(-1/3)), recombines the
 * algebraic part (Cancel[Together] with Numerator/Denominator Expand, when
 * it strictly shrinks), and distributes / sign-normalises inverse-trig
 * arguments.  Non-radical and partially-closed (nested Integrate) results
 * are returned unchanged.  TAKES OWNERSHIP of `r` and returns it (or a
 * derived tree). */
Expr* intsimp_finalize(Expr* r, Expr* x);

#endif /* MATHILDA_INTSIMP_H */
