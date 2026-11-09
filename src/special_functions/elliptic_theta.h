/* Mathilda -- EllipticTheta / EllipticThetaPrime (nome q). See elliptic_theta.c. */
#ifndef MATHILDA_ELLIPTIC_THETA_H
#define MATHILDA_ELLIPTIC_THETA_H

#include <stdbool.h>
#include "expr.h"

/* theta_a(u, q) (prime = 0) or its u-derivative (prime = 1), a = 1..4, real or
 * complex parts: the helper behind the scalar builtins, the NDArray / packed
 * map and the Compile[] VM. Declines only outside |q| < 1 or the double range. */
bool sf_machine_theta(int a, int prime, double ur, double ui, double qr, double qi,
                      double* rr, double* ri);

void elliptic_theta_init(void);

#endif /* MATHILDA_ELLIPTIC_THETA_H */
