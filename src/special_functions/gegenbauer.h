/* Mathilda -- GegenbauerC, the Gegenbauer (ultraspherical) polynomials and
 * functions C_n^lambda(z). See gegenbauer.c. */
#ifndef MATHILDA_GEGENBAUER_H
#define MATHILDA_GEGENBAUER_H

#include <stdbool.h>
#include "expr.h"

/* The machine kernels, shared by the scalar builtin, the NDArray / packed map
 * and the Compile[] VM (bit for bit). Real n, lambda, z (sf_machine_gegenbauer);
 * any complex part (sf_machine_gegenbauer_c). Decline only for a value outside
 * the normal double range (arb answers every point the double algorithms
 * cannot certify, inside the kernel). */
bool sf_machine_gegenbauer(double n, double lam, double x, double* out);
bool sf_machine_gegenbauer_c(double nr, double ni, double lr, double li,
                             double xr, double xi, double* rr, double* ri);

/* The two-argument form GegenbauerC[n, z] = 2 ChebyshevT[n, z] / n. */
bool sf_machine_gegenbauer2_c(double nr, double ni, double xr, double xi,
                              double* rr, double* ri);

Expr* builtin_gegenbauer_c(Expr* res);
void gegenbauer_init(void);

#endif /* MATHILDA_GEGENBAUER_H */
