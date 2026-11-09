/* Mathilda -- ParabolicCylinderD, Whittaker's parabolic cylinder function D_nu(z).
 * See parabolic_cylinder.c for the definitions and the algorithms. */
#ifndef MATHILDA_PARABOLIC_CYLINDER_H
#define MATHILDA_PARABOLIC_CYLINDER_H

#include <stdbool.h>
#include "expr.h"

/* The machine kernels: real order nu, real z (sf_machine_pcfd) or complex z
 * (sf_machine_pcfd_c). Shared by the scalar builtin, the NDArray / packed
 * element-wise map and the Compile[] VM, so every surface agrees bit for bit.
 * Return false ("decline") when the double algorithms cannot certify ~machine
 * accuracy, or the value leaves the normal double range; the caller then takes
 * the arbitrary-precision (arb) path. */
bool sf_machine_pcfd(double nu, double z, double* out);
bool sf_machine_pcfd_c(double nu, double zr, double zi, double* rr, double* ri);
/* Complex order as well (nu = nr + i ni): the arb machine twin when ni != 0. */
bool sf_machine_pcfd_cc(double nr, double ni, double zr, double zi, double* rr, double* ri);

Expr* builtin_parabolic_cylinder_d(Expr* res);

/* Register ParabolicCylinderD (attributes, docstring). Called from core_init(). */
void parabolic_cylinder_init(void);

#endif /* MATHILDA_PARABOLIC_CYLINDER_H */
