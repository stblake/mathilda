#ifndef MATHILDA_RIEMANNR_H
#define MATHILDA_RIEMANNR_H

#include <stdbool.h>
#include "expr.h"

/* RiemannR[x] -- Riemann's prime-counting function, by the Gram series
 * 1 + Sum_k Log[x]^k/(k k! Zeta[k+1]) (src/special_functions/riemannr.c). */
Expr* builtin_riemannr(Expr* res);
/* Internal`RiemannRPrime[z]: R'(z) at an inexact z (what Derivative[1][RiemannR]
 * evaluates to numerically). */
Expr* builtin_riemannr_prime(Expr* res);
void  riemannr_init(void);

/* The machine kernels shared by the scalar builtin, the NDArray / packed map
 * and Compile[].  sf_riemannr_r: real x >= 1/e (declines below, where the
 * alternating series cancels; the MPFR path answers).  sf_riemannr_c: complex
 * z with |Log z| <= 3/2, declining positive reals.  riemannr_machine_complex:
 * the unary-kernel complex ABI (both arms). */
bool sf_riemannr_r(double x, double* out);
bool sf_riemannr_c(double zr, double zi, double* ore, double* oim);
bool riemannr_machine_complex(double are, double aim, double* ore, double* oim);

#endif /* MATHILDA_RIEMANNR_H */
