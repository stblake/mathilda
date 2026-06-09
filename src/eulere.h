#ifndef EULERE_H
#define EULERE_H

#include "expr.h"

/* EulerE[n]     -- the Euler number E_n (exact integer).
 * EulerE[n, x]  -- the Euler polynomial E_n(x).
 *
 * Attributes: Listable, Protected.
 *
 * Evaluation strategy (see eulere.c for detail):
 *   - EulerE[n] for a non-negative integer n reduces to the exact integer
 *     E_n via the standard recurrence (lazily cached with GMP). Odd n give 0,
 *     E_0 = 1, E_2 = -1, E_4 = 5, ...
 *   - An inexact integer-valued first argument (Real / MPFR) evaluates the
 *     exact value numerically at machine or arbitrary precision.
 *   - EulerE[n, x] for a non-negative integer n builds the degree-n polynomial
 *     in expanded monomial form with exact rational coefficients and evaluates
 *     it: symbolic x stays symbolic, inexact x evaluates numerically.
 *   - EulerE[n, 1/2] with symbolic n folds to 2^-n EulerE[n].
 *   - Everything else stays symbolic (return NULL). */
Expr* builtin_eulere(Expr* res);

void eulere_init(void);

#endif /* EULERE_H */
