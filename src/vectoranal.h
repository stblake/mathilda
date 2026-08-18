/* vectoranal.h -- Vector-analysis differential operators.
 *
 * Grad, Div, Curl, Laplacian: gradient, divergence, curl, and Laplacian of
 * scalar / vector / tensor fields, in Cartesian coordinates (2-arg form) and in
 * the orthogonal coordinate charts "Cartesian", "Polar", "Cylindrical",
 * "Spherical" (3-arg form). All are assembled on top of the derivative engine
 * D (src/calculus/deriv.c); see vectoranal.c for the mathematical design.
 */
#ifndef VECTORANAL_H
#define VECTORANAL_H

#include "expr.h"

Expr* builtin_grad(Expr* res);
Expr* builtin_div(Expr* res);
Expr* builtin_curl(Expr* res);
Expr* builtin_laplacian(Expr* res);

/* Register Grad, Div, Curl, Laplacian (Protected) in the symbol table. */
void vectoranal_init(void);

#endif /* VECTORANAL_H */
