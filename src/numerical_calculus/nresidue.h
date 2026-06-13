/*
 * nresidue.h — NResidue[expr, {z, z0}, opts]
 *
 * Numerically estimates the residue of `expr` at z = z0 (the coefficient of
 * (z - z0)^-1 in the Laurent expansion) by contour integration in the
 * complex plane. Backed by the periodic-trapezoidal core in quadrature.h.
 *
 * Attributes: Protected (list-threading over arg 1 is handled manually so
 * the {z, z0} spec is never split — see nresidue.c).
 */
#ifndef MATHILDA_NRESIDUE_H
#define MATHILDA_NRESIDUE_H

#include "expr.h"

Expr* builtin_nresidue(Expr* res);
void  nresidue_init(void);

#endif /* MATHILDA_NRESIDUE_H */
