#ifndef GRUNTZ_H
#define GRUNTZ_H

#include "expr.h"

/*
 * Gruntz's algorithm for computing limits of exp-log functions.
 *
 * This is a faithful C port of the most-rapidly-varying (mrv) expansion
 * method from Dominik Gruntz's 1996 ETH PhD thesis "On Computing Limits in a
 * Symbolic Manipulation System" (Chapter 3 + Appendix A), following the
 * proven structure of SymPy's gruntz.py.
 *
 * The whole function is expanded as a power series in its most rapidly
 * varying subexpression w -> 0+, which structurally avoids the intermediate
 * expression swell that bottom-up series approaches suffer from mutual
 * cancellation.
 *
 * gruntz_limit computes  lim_{x -> point, direction} f.
 *
 *   x        the limit variable (a Symbol)
 *   point    the limit point: Infinity, -Infinity, or a finite expression
 *   dir      +1 = from above (LIMIT_DIR_FROMABOVE), -1 = from below,
 *            0 = two-sided (only meaningful at a finite point).
 *   depth    incoming recursion depth from the caller (for the guard).
 *
 * Inputs are borrowed (not freed). Returns a freshly-allocated owned Expr*
 * on success, or NULL when the engine cannot decide (non-exp-log input, an
 * undecidable sign/zero oracle, or an internal give-up) -- the caller then
 * leaves Limit unevaluated, per Mathilda's convention.
 */
Expr* gruntz_limit(Expr* f, Expr* x, Expr* point, int dir, int depth);

/* No symbol-table registration is required: the engine is reached through
 * limit.c's cascade. Provided for symmetry with the other *_init hooks. */
void gruntz_init(void);

#endif /* GRUNTZ_H */
