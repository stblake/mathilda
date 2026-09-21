/* algebraicintegerq.h — AlgebraicIntegerQ[x] predicate.
 *
 * AlgebraicIntegerQ[x] yields True if x is an algebraic integer (a root of a
 * monic polynomial with integer coefficients) and False otherwise.  Decided
 * exactly via FLINT's qqbar: x is an algebraic integer iff its primitive integer
 * minimal polynomial is monic.  Rational integers are algebraic integers;
 * non-integer rationals and non-algebraic quantities (Pi, a free symbol, a list)
 * are not.
 *
 * Attributes: Protected.  (Not Listable: a list is not an algebraic integer, so
 * AlgebraicIntegerQ[{..}] is False; map with /@ to test elements.)  Requires
 * FLINT — in a USE_FLINT=0 build the head is left unevaluated.
 *
 * Ownership: standard contract — a fresh True/False Expr on a decision, or NULL
 * to leave res unevaluated (bad arity, or FLINT compiled out).
 */
#ifndef MATHILDA_ALGEBRAICINTEGERQ_H
#define MATHILDA_ALGEBRAICINTEGERQ_H

#include "expr.h"

Expr* builtin_algebraicintegerq(Expr* res);
void algebraicintegerq_init(void);

#endif /* MATHILDA_ALGEBRAICINTEGERQ_H */
