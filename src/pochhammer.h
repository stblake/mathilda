#ifndef POCHHAMMER_H
#define POCHHAMMER_H

#include "expr.h"

/* Pochhammer[a, n] -- the Pochhammer symbol (rising factorial)
 *
 *   Pochhammer[a, n] = a (a+1) ... (a+n-1) = Gamma[a+n] / Gamma[a]
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * Evaluation strategy (see pochhammer.c for detail):
 *   - Exact integer n (|n| <= cap) expands to a Times product of n linear
 *     factors, which collapses to an exact BigInt / Rational for numeric a,
 *     preserves MPFR precision and complex arithmetic, and stays a symbolic
 *     polynomial product for symbolic a.
 *   - Otherwise, for numeric a and n, evaluates Gamma[a+n]/Gamma[a], reusing
 *     the Gamma builtin's exact, machine, MPFR and complex paths.
 *   - Everything else stays symbolic. */
Expr* builtin_pochhammer(Expr* res);

void pochhammer_init(void);

#endif /* POCHHAMMER_H */
