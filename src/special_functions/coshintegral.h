#ifndef COSHINTEGRAL_H
#define COSHINTEGRAL_H

#include "expr.h"

/* CoshIntegral[z] -- the hyperbolic cosine integral
 * Chi(z) = EulerGamma + Log[z] + Int_0^z (Cosh[t] - 1)/t dt.
 * Has a logarithmic singularity at 0 and a branch cut along the negative real
 * axis (from -Infinity to 0). */
Expr* builtin_coshintegral(Expr* res);
void  coshintegral_init(void);

#endif /* COSHINTEGRAL_H */
