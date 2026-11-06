/* Mathilda -- RiemannSiegelZ and RiemannSiegelTheta.
 *
 *   RiemannSiegelTheta[t]   theta(t) = Im LogGamma(1/4 + i t/2) - t Log(Pi)/2
 *   RiemannSiegelZ[t]       Z(t) = Exp(i theta(t)) Zeta(1/2 + i t), real for
 *                           real t, with |Z(t)| = |Zeta(1/2 + i t)|
 *
 * See riemann_siegel.c for the evaluation strategy.
 */
#ifndef MATHILDA_RIEMANN_SIEGEL_H
#define MATHILDA_RIEMANN_SIEGEL_H

#include <stdbool.h>
#include "expr.h"

Expr* builtin_riemann_siegel_z(Expr* res);
Expr* builtin_riemann_siegel_theta(Expr* res);

/* The machine helpers behind the scalar builtins, the NDArray kernels and the
 * Compile[] lowering (one helper per function, so the surfaces agree bit for
 * bit). Real and complex t; false = declined (non-finite, FLINT absent). */
bool sf_machine_riemann_siegel_z(double t, double* out);
bool sf_machine_riemann_siegel_z_c(double tr, double ti, double* rr, double* ri);
bool sf_machine_riemann_siegel_theta(double t, double* out);
bool sf_machine_riemann_siegel_theta_c(double tr, double ti, double* rr, double* ri);

/* Register both builtins, their attributes and docstrings. */
void riemann_siegel_init(void);

#endif /* MATHILDA_RIEMANN_SIEGEL_H */
