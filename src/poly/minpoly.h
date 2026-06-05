/*
 * minpoly.h
 * ---------
 * MinimalPolynomial[s, x] — the minimal (defining) polynomial over Q of an
 * algebraic number s, returned as a primitive integer-coefficient polynomial
 * in x with positive leading coefficient.
 *
 *   MinimalPolynomial[s, x]                lowest-degree integer polynomial in
 *                                          x having the algebraic number s as a
 *                                          root (content 1, positive lead).
 *   MinimalPolynomial[s]                   the same, as a pure function
 *                                          (e.g. 1 - 10 #1^2 + #1^4 &).
 *   MinimalPolynomial[s, x, Extension->a]  the characteristic polynomial of
 *                                          s in Q(a) over Q(a) — m_s(x)^(d/e)
 *                                          where d = [Q(a):Q], e = [Q(s):Q],
 *                                          provided s lies in Q(a).
 *
 * Threads over lists (Listable).  The algebraic number may be built from
 * integers/rationals, radicals (Sqrt / Power[_, Rational]), the imaginary
 * unit, roots of unity (Power[E, I Pi r]), Root[] objects, and the field
 * operations Plus/Times/Power.
 *
 * Algorithm: each algebraic atom is replaced by a fresh auxiliary variable
 * with a polynomial defining relation; the auxiliary variables are eliminated
 * from (x - value) by repeated Resultant; the result is made primitive and
 * factored over Z, and the unique irreducible factor vanishing at s (selected
 * by high-precision numeric evaluation) is returned.
 */
#ifndef MINPOLY_H
#define MINPOLY_H

#include "expr.h"

/* Builtin entry point.  Takes ownership of res per the standard contract:
 * returns a fresh Expr* on success, or NULL (leaving the call unevaluated)
 * for non-algebraic input or a malformed call. */
Expr* builtin_minimalpolynomial(Expr* res);

/* Register MinimalPolynomial with the symbol table and set its attributes. */
void minpoly_init(void);

#endif /* MINPOLY_H */
