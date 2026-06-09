#ifndef GAMMA_H
#define GAMMA_H

#include "expr.h"

/* Gamma[z]          -- Euler gamma function  Gamma(z)
 * Gamma[a, z]       -- upper incomplete gamma Gamma(a, z) = Int_z^Inf t^(a-1) e^-t dt
 * Gamma[a, z0, z1]  -- generalized incomplete = Gamma(a, z0) - Gamma(a, z1)
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * Evaluation strategy (see gamma.c for detail):
 *   - Exact integer / half-integer arguments reduce to (z-1)! via the
 *     existing Factorial machinery (exact / BigInt / rational*Sqrt[Pi]).
 *   - Machine-precision real -> libm tgamma; complex -> Lanczos.
 *   - Arbitrary-precision real -> MPFR mpfr_gamma (and mpfr_gamma_inc for
 *     the incomplete form).
 *   - Everything else stays symbolic. */
Expr* builtin_gamma(Expr* res);

void gamma_init(void);

#endif /* GAMMA_H */
