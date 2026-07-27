#ifndef COSINTEGRAL_H
#define COSINTEGRAL_H

#include "expr.h"
#include <stdbool.h>

/* CosIntegral[z] -- the cosine integral  Ci(z) = -Int_z^Inf Cos[t]/t dt.
 * Has a logarithmic singularity at 0 and a branch cut along the negative
 * real axis (from -Infinity to 0). */
Expr* builtin_cosintegral(Expr* res);
void  cosintegral_init(void);


/* CosIntegral for a machine complex argument, in the shared kernel ABI.  Returns false
 * to DECLINE — at the singularity, and wherever the ascending series has lost
 * too much of the answer to cancellation for a double to carry it (the MPFR
 * path answers there). */
bool cosintegral_machine_complex(double are, double aim, double* ore, double* oim);

#endif /* COSINTEGRAL_H */
