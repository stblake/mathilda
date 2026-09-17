#ifndef FACPOLY_LIST_H
#define FACPOLY_LIST_H

#include "expr.h"

/* FactorList[poly]
 * FactorList[poly, GaussianIntegers -> True]
 * FactorList[poly, Extension -> {a1, a2, ...}]
 *
 * Gives a list of the irreducible factors of `poly` together with their
 * exponents, as {factor, exponent} pairs.  A thin wrapper over Factor: it
 * factors via `Factor[poly, opts...]` (options are forwarded verbatim) and
 * splits the product into pairs.
 *
 * The first element is always the overall numerical factor {c, 1} (it is
 * {1, 1} when there is no numerical factor).  Denominator factors of a
 * rational function appear with negative exponents.
 *
 * Wrong arg count emits `FactorList::argx`; extra non-option arguments beyond
 * position 1 emit `FactorList::nonopt`; both return NULL (call unevaluated). */
Expr* builtin_factorlist(Expr* res);

/* FactorSquareFreeList[poly]
 * FactorSquareFreeList[poly, Extension -> ext]
 *
 * Gives a list of the square-free factors of `poly` together with their
 * multiplicities, as {factor, exponent} pairs.  A thin wrapper over
 * FactorSquareFree: it forwards all arguments (poly + the Extension option)
 * verbatim and splits the resulting product into pairs.  The first element is
 * the overall numerical factor {c, 1} (it is {1, 1} when there is none).
 * FactorSquareFreeList[] (no argument) returns NULL (call unevaluated). */
Expr* builtin_factorsquarefreelist(Expr* res);

void factorlist_init(void);

#endif /* FACPOLY_LIST_H */
