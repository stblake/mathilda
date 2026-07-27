#ifndef SINHINTEGRAL_H
#define SINHINTEGRAL_H

#include "expr.h"
#include <stdbool.h>

/* SinhIntegral[z] -- the hyperbolic sine integral  Shi(z) = Int_0^z Sinh[t]/t dt.
 * An entire, odd function with no branch cuts. */
Expr* builtin_sinhintegral(Expr* res);
void  sinhintegral_init(void);


/* SinhIntegral for a machine complex argument, in the shared kernel ABI.  Returns false
 * to DECLINE — at the singularity, and wherever the ascending series has lost
 * too much of the answer to cancellation for a double to carry it (the MPFR
 * path answers there). */
bool sinhintegral_machine_complex(double are, double aim, double* ore, double* oim);

#endif /* SINHINTEGRAL_H */
