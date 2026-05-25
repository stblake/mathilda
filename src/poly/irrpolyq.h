/*
 * irrpolyq.h -- IrreduciblePolynomialQ[poly, opts].
 *
 *   IrreduciblePolynomialQ[poly]
 *     True iff poly is an irreducible polynomial over Q.
 *
 *   IrreduciblePolynomialQ[poly, Extension -> α | {α_1, ..., α_n}]
 *     Tests irreducibility over Q(α) (or Q(α_1, ..., α_n)).
 *
 *   IrreduciblePolynomialQ[poly, Extension -> Automatic]
 *     Tests irreducibility over Q extended by every algebraic-number
 *     coefficient appearing in poly.
 *
 *   IrreduciblePolynomialQ[poly, Extension -> All]
 *     Tests absolute irreducibility over the complex numbers.
 *
 *   IrreduciblePolynomialQ[poly, GaussianIntegers -> True]
 *     Tests irreducibility over the Gaussian rationals Q(i).
 *
 * Attributes: Listable, Protected.
 * Returns True / False on a structurally valid call; emits ::argx or
 * ::nonopt on arg-shape errors and returns NULL (leaving the call
 * unevaluated, matching Mathematica's surface behaviour).
 */

#ifndef MATHILDA_POLY_IRRPOLYQ_H
#define MATHILDA_POLY_IRRPOLYQ_H

#include "expr.h"

Expr* builtin_irreduciblepolynomialq(Expr* res);
void  irrpolyq_init(void);

#endif /* MATHILDA_POLY_IRRPOLYQ_H */
