/* algebraicnumbertrace.h — AlgebraicNumberTrace[a].
 *
 * AlgebraicNumberTrace[a] gives the field trace of the algebraic number a: the
 * sum of the conjugates of a over the rationals, i.e. the sum of the roots of
 * a's minimal polynomial.  It accepts the same surface as the other
 * AlgebraicNumber* heads: integers, rationals, radicals, the named constant
 * GoldenRatio, Root objects, and AlgebraicNumber objects.
 *
 * AlgebraicNumberTrace[a, Extension -> theta] gives the relative trace over the
 * field Q(theta); a must lie in Q(theta).  The value is computed exactly from
 * minimal polynomials (flint_qqbar_algebraic_number_trace): the absolute trace is
 * -(coeff of x^{deg-1}) / (leading coeff), and the relative trace is the absolute
 * trace scaled by the tower index [Q(theta):Q(a)].
 *
 * Attributes: Listable, Protected.  (Listable, so it threads element-wise over a
 * list of algebraic numbers, with any trailing Extension option repeated.)
 * Requires FLINT — in a USE_FLINT=0 build the head is left unevaluated.
 *
 * Ownership: standard contract — a fresh Integer/Rational on success, or NULL to
 * leave res unevaluated (bad arity, a non-algebraic argument, an argument not in
 * Q(theta), or FLINT compiled out).  A non-algebraic argument emits
 * AlgebraicNumberTrace::nalg; an argument not in the extension emits
 * AlgebraicNumberTrace::ext.
 */
#ifndef MATHILDA_ALGEBRAICNUMBERTRACE_H
#define MATHILDA_ALGEBRAICNUMBERTRACE_H

#include "expr.h"

Expr* builtin_algebraicnumbertrace(Expr* res);
void algebraicnumbertrace_init(void);

#endif /* MATHILDA_ALGEBRAICNUMBERTRACE_H */
