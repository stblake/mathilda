#ifndef SINHINTEGRAL_H
#define SINHINTEGRAL_H

#include "expr.h"

/* SinhIntegral[z] -- the hyperbolic sine integral  Shi(z) = Int_0^z Sinh[t]/t dt.
 * An entire, odd function with no branch cuts. */
Expr* builtin_sinhintegral(Expr* res);
void  sinhintegral_init(void);

#endif /* SINHINTEGRAL_H */
