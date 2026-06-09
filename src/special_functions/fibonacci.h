/* Mathilda -- Fibonacci numbers and Fibonacci polynomials.
 *
 *   Fibonacci[n]     -- the nth Fibonacci number F_n.
 *   Fibonacci[n, x]  -- the Fibonacci polynomial F_n(x), the coefficient of
 *                       t^n in the expansion of t / (1 - x t - t^2).
 *
 * Exact integer orders are computed with GMP (fast doubling); non-integer
 * or inexact orders use the generalized closed form and the existing numeric
 * machinery (machine double / MPFR). See fibonacci.c for details.
 */
#ifndef FIBONACCI_H
#define FIBONACCI_H

#include "expr.h"

/* Fibonacci[n] / Fibonacci[n, x]. Returns a freshly allocated Expr on
 * success, or NULL to leave the call unevaluated (symbolic order). */
Expr* builtin_fibonacci(Expr* res);

/* Register Fibonacci with the symbol table. */
void fibonacci_init(void);

#endif /* FIBONACCI_H */
