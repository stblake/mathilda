#ifndef FRESNEL_H
#define FRESNEL_H

#include "expr.h"

/* Fresnel integrals (Pi/2-normalized):
 *   FresnelC[z] = Integral_0^z Cos[Pi t^2 / 2] dt
 *   FresnelS[z] = Integral_0^z Sin[Pi t^2 / 2] dt
 * Both are entire, odd functions with no branch cuts. */
Expr* builtin_fresnelc(Expr* res);
Expr* builtin_fresnels(Expr* res);

/* Registers FresnelC and FresnelS (both builtins) with attributes and
 * interns their symbols. Docstrings live in info.c. */
void fresnel_init(void);

#endif /* FRESNEL_H */
