/* algebraicnumberdenominator.h — AlgebraicNumberDenominator[a].
 *
 * AlgebraicNumberDenominator[a] gives the smallest positive integer n such that
 * n a is an algebraic integer.  It accepts the same surface as AlgebraicIntegerQ:
 * rationals, radicals (1/Sqrt[3], (1+3I)^(-1/3), ...), Root objects, and
 * AlgebraicNumber objects.  The value is computed exactly from the minimal
 * polynomial via a per-prime valuation (flint_qqbar_algebraic_number_denominator)
 * -- not the raw min-poly leading coefficient, which over-counts.
 *
 * Attributes: Listable, Protected.  (Listable, so it threads element-wise over a
 * list of algebraic numbers.)  Requires FLINT — in a USE_FLINT=0 build the head is
 * left unevaluated.
 *
 * Ownership: standard contract — a fresh Integer/BigInt on success, or NULL to
 * leave res unevaluated (bad arity, a non-algebraic argument, or FLINT compiled
 * out).  A non-algebraic argument also emits AlgebraicNumberDenominator::nalg.
 */
#ifndef MATHILDA_ALGEBRAICNUMBERDENOMINATOR_H
#define MATHILDA_ALGEBRAICNUMBERDENOMINATOR_H

#include "expr.h"

Expr* builtin_algebraicnumberdenominator(Expr* res);
void algebraicnumberdenominator_init(void);

#endif /* MATHILDA_ALGEBRAICNUMBERDENOMINATOR_H */
