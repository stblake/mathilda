/* symrad.h -- Together / Cancel over Q(variables)(Sqrt[b]) for a polynomial
 * radicand b.
 *
 * A square root of a polynomial in the variables (the Hessian roots
 * (1 + k +- Sqrt[-3 + 6 k - 3 k^2])/(2 - 2 k + 2 k^2) of a cubic with a
 * parameter k) is algebraic over the rational function field Q(vars), not over
 * Q: the field the expression lives in is K = Q(vars)[s]/(s^2 - b).  No other
 * engine of Together / Cancel covers it -- the number-field engines need
 * constant algebraic numbers, the parametric engine a bare-symbol radicand
 * (Q(k)(Sqrt[k]) is the rational function field Q(Sqrt[k])) -- so such input fell
 * to the multivariate gcds of the classical cascade over Q[vars, Sqrt[b]], which
 * do not see s^2 = b at all and swell without bound (seconds to minutes for a
 * single Cancel of a few factors).
 *
 * K is a field whenever the squarefree part of b is not a square, and then
 * every element has the canonical form (A + B Sqrt[b])/C with A, B, C
 * polynomials free of Sqrt[b]: combine with s a fresh symbol (one Together over
 * Q(vars, s)), reduce numerator and denominator modulo s^2 - b, multiply by the
 * conjugate of the denominator, so that only the norm d0^2 - d1^2 b, a
 * polynomial in the variables, remains below, and reduce each coordinate in
 * lowest terms.  Everything runs over Q(vars) (or over Q(constants)(vars) when
 * the input also carries constant algebraic numbers), where the factored
 * rational-function engine is fast. */
#ifndef SYMRAD_H
#define SYMRAD_H

#include "expr.h"

/* `head` is SYM_Together or SYM_Cancel; `e` the expression (the Extension ->
 * Automatic option already stripped).  Returns the canonical form over
 * Q(vars)(Sqrt[b]) when `e` carries square roots of exactly one non-constant
 * polynomial radicand b (not a bare symbol, which the parametric engine owns),
 * the squarefree part of b is not a square, and nothing else beyond rational
 * numbers, variables and constant algebraic numbers -- or NULL to decline.  A
 * radicand whose squarefree part is a constant c (-3 (k - 1)^2) is declined
 * when Sqrt[c] lies in the field of the input's constant algebraic numbers
 * (Sqrt[-3] in Q((-1)^(1/3))): Sqrt[b]/(k - 1) then squares to a square there,
 * and K is a ring with zero divisors.  Never mutates `e`. */
Expr* symrad_normalize(const Expr* e, const char* head);

#endif /* SYMRAD_H */
