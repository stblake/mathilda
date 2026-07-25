/* Mathilda — compile a linear method-of-lines RHS into a numeric operator.
 *
 * If every reduced right-hand side f_i(t, Y) is affine in the reduced-state
 * symbols with constant coefficients — i.e. dY/dt = A·Y + s(t) — this builds the
 * constant matrix A (with its bandwidth) and the forcing s(t), so the driver can
 * evaluate the RHS by a matrix–vector product and use A as the exact Jacobian.
 * Returns NULL when the system is nonlinear or the coefficients depend on t
 * (in which case the symbolic sampler is used instead). */
#ifndef NDSOLVE_OPERATOR_H
#define NDSOLVE_OPERATOR_H

#include "ndsolve_common.h"

NdOperator* nd_operator_try_build(NdProblem* P);

#endif /* NDSOLVE_OPERATOR_H */
