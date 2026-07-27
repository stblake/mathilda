#ifndef LOGINTEGRAL_H
#define LOGINTEGRAL_H

#include "expr.h"
#include <stdbool.h>

/* LogIntegral[z] -- the logarithmic integral li(z) = PV Int_0^z dt/ln t,
 * with a branch cut along (-Infinity, +1). Computed as Ei(Log z). */
Expr* builtin_logintegral(Expr* res);
void  logintegral_init(void);


/* li(z) for a machine complex argument, in the shared kernel ABI.  Returns
 * false to DECLINE (at z = 0 and z = 1, and wherever Ei declines). */
bool logintegral_machine_complex(double are, double aim, double* ore, double* oim);

#endif /* LOGINTEGRAL_H */
