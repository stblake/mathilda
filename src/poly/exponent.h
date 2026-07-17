#ifndef EXPONENT_H
#define EXPONENT_H

#include "expr.h"

/* Exponent[expr, form]        -- maximum power of `form` in the expanded expr
 * Exponent[expr, form, h]     -- applies h to the set of exponents (default Max)
 *
 * `form` may be a symbol, a kernel (e.g. Sin[x], Sqrt[2]) or a product of
 * terms.  Exponent is purely syntactic: it expands `expr` but does not attempt
 * to recognise zero coefficients.  Rational and symbolic exponents are
 * supported (e.g. Exponent[x^(n+1) + 2 Sqrt[x] + 1, x] = Max[1/2, 1 + n]).
 *
 *   Exponent[0, x]            = -Infinity   (empty exponent set, h = Max)
 *   Exponent[expr, x, List]   = sorted, de-duplicated list of exponents
 *   Exponent[expr, x, Min]    = lowest exponent
 *
 * The Listable attribute threads Exponent[expr, {f1, f2, ...}] into the list of
 * per-form exponents.  Wrong arg count emits `Exponent::argt`; extra non-option
 * arguments beyond position 3 emit `Exponent::nonopt`; both return NULL. */
Expr* builtin_exponent(Expr* res);

void exponent_init(void);

#endif /* EXPONENT_H */
