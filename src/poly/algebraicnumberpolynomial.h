/* algebraicnumberpolynomial.h — AlgebraicNumberPolynomial[a, x] builtin.
 *
 * Gives the polynomial in x that defines the AlgebraicNumber object a: for
 * a = AlgebraicNumber[theta, {c0, c1, ..., cn}], returns c0 + c1 x + ... +
 * cn x^n, so that a is recovered by the substitution x -> theta. An integer or
 * rational a is the constant polynomial and is returned unchanged.
 *
 * The coefficient vector is already stored in the AlgebraicNumber object, so
 * this is a purely structural read + polynomial build — no FLINT is required
 * (and none is used): the head works even in a USE_FLINT=0 build (though such a
 * build never produces AlgebraicNumber objects to feed it).
 *
 * Ownership: builtin_algebraicnumberpolynomial follows the standard contract —
 * returns a fresh Expr on success (the evaluator frees res) or NULL to leave
 * res as is (a non-AlgebraicNumber, non-rational argument, after emitting the
 * AlgebraicNumberPolynomial::naobj message).
 */
#ifndef MATHILDA_ALGEBRAICNUMBERPOLYNOMIAL_H
#define MATHILDA_ALGEBRAICNUMBERPOLYNOMIAL_H

#include "expr.h"

Expr* builtin_algebraicnumberpolynomial(Expr* res);
void algebraicnumberpolynomial_init(void);

#endif /* MATHILDA_ALGEBRAICNUMBERPOLYNOMIAL_H */
