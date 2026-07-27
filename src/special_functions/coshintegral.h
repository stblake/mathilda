#ifndef COSHINTEGRAL_H
#define COSHINTEGRAL_H

#include "expr.h"
#include <stdbool.h>

/* CoshIntegral[z] -- the hyperbolic cosine integral
 * Chi(z) = EulerGamma + Log[z] + Int_0^z (Cosh[t] - 1)/t dt.
 * Has a logarithmic singularity at 0 and a branch cut along the negative real
 * axis (from -Infinity to 0). */
Expr* builtin_coshintegral(Expr* res);
void  coshintegral_init(void);


/* CoshIntegral for a machine complex argument, in the shared kernel ABI.  Returns false
 * to DECLINE — at the singularity, and wherever the ascending series has lost
 * too much of the answer to cancellation for a double to carry it (the MPFR
 * path answers there). */
bool coshintegral_machine_complex(double are, double aim, double* ore, double* oim);

#endif /* COSHINTEGRAL_H */
