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

/* Entry point bound to the Eliminate symbol in eliminate_init(). */
Expr* builtin_eliminate(Expr* res);

/* Register the builtin and set ATTR_PROTECTED. */
void eliminate_init(void);

#ifdef __cplusplus
}
#endif

#endif /* POLY_ELIMINATE_H */
