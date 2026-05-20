/*
 * limit.h -- Symbolic limit evaluation (Limit[f, x -> a], Direction, ...).
 *
 * This module provides a native C implementation of Mathematica's Limit
 * built-in, following the pipeline laid out in plans/limit_candidate_spec.md.
 * The public surface is deliberately small -- a single registration hook
 * (limit_init) plus the evaluator entry point (builtin_limit). Everything
 * else is internal; the full dispatch chain is documented in limit.c.
 *
 * Calling forms (per Mathematica):
 *
 *     Limit[f, x -> a]
 *     Limit[f, {x1 -> a1, ..., xn -> an}]           (iterated, rightmost first)
 *     Limit[f, {x1, ..., xn} -> {a1, ..., an}]      (multivariate)
 *     Limit[f, x -> a, Direction -> spec]
 *
 * Direction settings accepted by the front end:
 *   Reals, "TwoSided"         -- two-sided limit (default on the real line)
 *   "FromAbove" or -1         -- approach from above (Mathematica convention)
 *   "FromBelow" or +1         -- approach from below (Mathematica convention)
 *   Complexes                 -- limit taken over all complex directions
 *   Exp[I theta]              -- limit along the given complex ray
 *   {dir1, ..., dirn}         -- per-variable directions in the multivariate
 *                                form
 *
 * Note: Mathematica's Direction sign convention is the opposite of
 * MACSYMA's. The dispatcher in limit.c performs the sign flip before
 * dispatching into the internal ZERO_PLUS / ZERO_MINUS tags.
 *
 * Return values (summarised):
 *   finite expression / Infinity / -Infinity   -- limit exists
 *   Indeterminate                               -- provably no limit
 *   Interval[{lo, hi}]                         -- bounded, but no single value
 *   the original Limit[...] expression          -- could not determine
 */

#ifndef MATHILDA_LIMIT_H
#define MATHILDA_LIMIT_H

#include "expr.h"

/* Register Limit with the symbol table. Idempotent. */
void limit_init(void);

/* Built-in entry point. Follows the standard Mathilda contract:
 *   - caller owns `res`;
 *   - return NULL to leave `res` unevaluated (caller retains ownership);
 *   - return a freshly allocated Expr* to replace `res` (caller frees `res`).
 */
Expr* builtin_limit(Expr* res);

#endif /* MATHILDA_LIMIT_H */
