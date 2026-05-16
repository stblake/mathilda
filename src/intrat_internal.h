/* intrat_internal.h — small private surface shared between intrat.c
 * and intsimp.c.  Not part of the public Integrate` package API; do
 * not include outside the rational-function integrator.
 *
 * Each helper returns a freshly-allocated Expr* the caller owns.  None
 * of them free their inputs — every helper makes its own copies.
 */

#ifndef MATHILDA_INTRAT_INTERNAL_H
#define MATHILDA_INTRAT_INTERNAL_H

#include <stdbool.h>
#include "expr.h"

/* canonic[expr] = Cancel[Together[expr]]. */
Expr* intrat_canonic(Expr* e);

/* Denominator[expr]. */
Expr* intrat_denominator(Expr* e);

/* TrueQ[FreeQ[expr, var]]. */
bool  intrat_freeq_test(Expr* expr, Expr* var);

#endif /* MATHILDA_INTRAT_INTERNAL_H */
