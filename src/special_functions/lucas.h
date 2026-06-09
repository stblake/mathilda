/* Mathilda -- Lucas numbers and Lucas polynomials.
 *
 *   LucasL[n]     -- the nth Lucas number L_n.
 *   LucasL[n, x]  -- the Lucas polynomial L_n(x), the coefficient of t^n in
 *                    the expansion of (2 - t x) / (1 - x t - t^2).
 *
 * Exact integer orders are computed with GMP (fast doubling, via the
 * identity L_m = 2 F_{m+1} - F_m); non-integer or inexact orders use the
 * generalized closed form and the existing numeric machinery (machine
 * double / MPFR). See lucas.c for details.
 */
#ifndef LUCAS_H
#define LUCAS_H

#include "expr.h"

/* LucasL[n] / LucasL[n, x]. Returns a freshly allocated Expr on success,
 * or NULL to leave the call unevaluated (symbolic order). */
Expr* builtin_lucasl(Expr* res);

/* Register LucasL with the symbol table. */
void lucas_init(void);

#endif /* LUCAS_H */
