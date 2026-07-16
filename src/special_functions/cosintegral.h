#ifndef COSINTEGRAL_H
#define COSINTEGRAL_H

#include "expr.h"

/* CosIntegral[z] -- the cosine integral  Ci(z) = -Int_z^Inf Cos[t]/t dt.
 * Has a logarithmic singularity at 0 and a branch cut along the negative
 * real axis (from -Infinity to 0). */
Expr* builtin_cosintegral(Expr* res);
void  cosintegral_init(void);

#endif /* COSINTEGRAL_H */
