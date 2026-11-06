#ifndef BESSELJZERO_H
#define BESSELJZERO_H

#include "expr.h"
#include <stdbool.h>

/* BesselJZero[n, k] -- the k-th positive zero of BesselJ[n, x]. */
Expr* builtin_besseljzero(Expr* res);
void  besseljzero_init(void);

/* The machine kernel behind the scalar builtin, the NDArray / packed map and
 * Compile[]: the k-th positive zero of J_nu, nu > -1 (at most 10^4), k a
 * positive integer (at most 10^6) given as a double. False outside that. */
bool sf_bessel_j_zero_r(double nu, double k, double* out);
/* BesselYZero[n, k]: the same for Y_nu, nu >= 0. */
Expr* builtin_besselyzero(Expr* res);
bool sf_bessel_y_zero_r(double nu, double k, double* out);

#endif /* BESSELJZERO_H */
