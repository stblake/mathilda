#ifndef SURD_H
#define SURD_H

#include "expr.h"
#include <stdbool.h>

/* Surd[x, n] -- the real-valued n-th root of a real x: Sign[x] Abs[x]^(1/n)
 * for odd n, x^(1/n) for even n and x >= 0. */
Expr* builtin_surd(Expr* res);
void  surd_init(void);

/* The machine kernel behind the scalar builtin, the NDArray / packed map and
 * Compile[]: Surd at doubles. Declines (false) on a non-integer or zero n, an
 * even root of a negative x, 0 to a negative index, and overflow -- each of
 * which the exact path answers (with its message) instead. */
bool sf_surd_r(double x, double n, double* out);

/* CubeRoot[x] -- Surd[x, 3] under its own name; sf_cuberoot_r is its machine
 * kernel (scalar builtin, NDArray / packed map and Compile[] alike). */
Expr* builtin_cuberoot(Expr* res);
bool  sf_cuberoot_r(double x, double* out);

#endif /* SURD_H */
