#ifndef BERNOULLIB_H
#define BERNOULLIB_H

#include "expr.h"

/* BernoulliB[n]     -- the Bernoulli number B_n (exact rational).
 * BernoulliB[n, x]  -- the Bernoulli polynomial B_n(x).
 *
 * Attributes: Listable, Protected.
 *
 * Evaluation strategy (see bernoullib.c for detail):
 *   - BernoulliB[n] for a non-negative integer n reduces to the exact
 *     rational B_n via the standard recurrence (lazily cached with GMP).
 *     Odd n > 1 give 0, B_0 = 1, B_1 = -1/2.
 *   - An inexact integer-valued first argument (Real / MPFR) evaluates the
 *     exact value numerically at machine or arbitrary precision.
 *   - BernoulliB[n, x] for a non-negative integer n builds the degree-n
 *     polynomial Sum_j C(n,j) B_{n-j} x^j with exact rational coefficients
 *     and evaluates it: symbolic x stays symbolic, inexact x evaluates
 *     numerically.
 *   - Everything else stays symbolic (return NULL). */
Expr* builtin_bernoullib(Expr* res);

void bernoullib_init(void);

#endif /* BERNOULLIB_H */
