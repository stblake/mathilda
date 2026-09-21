/* algebraicnumbernorm.h — AlgebraicNumberNorm[a].
 *
 * AlgebraicNumberNorm[a] gives the field norm of the algebraic number a: the
 * product of the conjugates of a over the rationals, i.e. the product of the
 * roots of a's minimal polynomial.  It accepts the same surface as the other
 * AlgebraicNumber* heads: integers, rationals, radicals, the named constant
 * GoldenRatio, Root objects, and AlgebraicNumber objects.
 *
 * AlgebraicNumberNorm[a, Extension -> theta] gives the relative norm over the
 * field Q(theta); a must lie in Q(theta).  The value is computed exactly from
 * minimal polynomials (flint_qqbar_algebraic_number_norm): the absolute norm is
 * (-1)^deg times the monic-over-Q constant term, and the relative norm is the
 * absolute norm raised to the tower index [Q(theta):Q(a)].
 *
 * Attributes: Listable, Protected.  (Listable, so it threads element-wise over a
 * list of algebraic numbers, with any trailing Extension option repeated.)
 * Requires FLINT — in a USE_FLINT=0 build the head is left unevaluated.
 *
 * Ownership: standard contract — a fresh Integer/Rational on success, or NULL to
 * leave res unevaluated (bad arity, a non-algebraic argument, an argument not in
 * Q(theta), or FLINT compiled out).  A non-algebraic argument emits
 * AlgebraicNumberNorm::nalg; an argument not in the extension emits
 * AlgebraicNumberNorm::ext.
 */
#ifndef MATHILDA_ALGEBRAICNUMBERNORM_H
#define MATHILDA_ALGEBRAICNUMBERNORM_H

#include "expr.h"

Expr* builtin_algebraicnumbernorm(Expr* res);
void algebraicnumbernorm_init(void);

#endif /* MATHILDA_ALGEBRAICNUMBERNORM_H */
