#ifndef GAMMA_H
#define GAMMA_H

#include "expr.h"
#include <stdbool.h>

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

/* Machine-complex Gamma, in the shared kernel ABI (returns false to decline, at
 * the poles and wherever the Lanczos series does not produce a finite value).
 *
 * Exposed so the Compile[] engine and the NDArray element-wise path call THE
 * SAME Lanczos series the interpreter does, rather than carrying a second copy:
 * a second copy is how one gets a branch-cut fix and the other does not. */
bool gamma_machine_complex(double are, double aim, double* ore, double* oim);

void gamma_init(void);

#endif /* GAMMA_H */
