/* Mathilda -- BetaRegularized[z, a, b] and BetaRegularized[z0, z1, a, b], the
 * regularized incomplete beta function. See betaregularized.c. */
#ifndef MATHILDA_BETAREGULARIZED_H
#define MATHILDA_BETAREGULARIZED_H

#include <stdbool.h>
#include "expr.h"

/* I_z(a, b) at real or complex (z, a, b): the helper behind the scalar builtin,
 * the NDArray / packed map and the Compile[] VM. The real case 0 <= z <= 1,
 * a, b > 0 runs the continued fraction in double; everything else, and a value
 * the double route cannot certify, is answered by arb at 53 bits. */
bool sf_machine_betareg(double zr, double zi, double ar, double ai, double br, double bi,
                        double* rr, double* ri);

/* I_z1(a, b) - I_z0(a, b), the four-argument form. */
bool sf_machine_betareg4(double z0r, double z0i, double z1r, double z1i, double ar, double ai,
                         double br, double bi, double* rr, double* ri);

void betaregularized_init(void);

#endif /* MATHILDA_BETAREGULARIZED_H */
