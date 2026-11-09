/* Mathilda -- HankelH1/H2 and the spherical Bessel family. See
 * hankel_spherical.c. */
#ifndef MATHILDA_HANKEL_SPHERICAL_H
#define MATHILDA_HANKEL_SPHERICAL_H

#include <stdbool.h>
#include "expr.h"

enum { HS_H1, HS_H2, HS_SJ, HS_SY, HS_SH1, HS_SH2, HS_COUNT };

/* The machine kernel for function `code` at (order, z), real or complex parts:
 * the helper behind the scalar builtins, the NDArray / packed map and the
 * Compile[] VM. Declines only for a pole or a value outside the double range. */
bool sf_machine_hankel_sph(int code, double nr, double ni, double zr, double zi,
                           double* rr, double* ri);

void hankel_spherical_init(void);

#endif /* MATHILDA_HANKEL_SPHERICAL_H */
