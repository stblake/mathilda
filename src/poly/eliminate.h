/* eliminate.h
 *
 * Eliminate[eqns, vars] -- user-facing variable elimination on systems of
 * polynomial equations.  Built as a thin driver over the existing lex-order
 * Buchberger engine in groebner.c (with GB_ORDER_ELIM for the elimination
 * block) and the principal-branch inverse-function rewriter described in
 * the changelog.
 */

#ifndef POLY_ELIMINATE_H
#define POLY_ELIMINATE_H

#include "expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/* When non-zero, Eliminate suppresses its advisory diagnostics
 * (::ifun, ::alg, ::nlin, ...).  Internal callers that drive Eliminate as a
 * private sub-step -- e.g. the DerivativeDivides integrator -- raise this
 * around their call so the user does not see branch-caveat messages for an
 * elimination they never asked for.  Save/restore it (a plain global, not
 * re-entrant). */
extern int eliminate_suppress_messages;

/* Entry point bound to the Eliminate symbol in eliminate_init(). */
Expr* builtin_eliminate(Expr* res);

/* Register the builtin and set ATTR_PROTECTED. */
void eliminate_init(void);

#ifdef __cplusplus
}
#endif

#endif /* POLY_ELIMINATE_H */
