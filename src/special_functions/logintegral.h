#ifndef LOGINTEGRAL_H
#define LOGINTEGRAL_H

#include "expr.h"

/* LogIntegral[z] -- the logarithmic integral li(z) = PV Int_0^z dt/ln t,
 * with a branch cut along (-Infinity, +1). Computed as Ei(Log z). */
Expr* builtin_logintegral(Expr* res);
void  logintegral_init(void);

#endif /* LOGINTEGRAL_H */
