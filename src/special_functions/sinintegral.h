#ifndef SININTEGRAL_H
#define SININTEGRAL_H

#include "expr.h"
#include <stdbool.h>

/* SinIntegral[z] -- the sine integral  Si(z) = Int_0^z Sin[t]/t dt.
 * An entire, odd function with no branch cuts. */
Expr* builtin_sinintegral(Expr* res);
void  sinintegral_init(void);


/* SinIntegral for a machine complex argument, in the shared kernel ABI.  Returns false
 * to DECLINE — at the singularity, and wherever the ascending series has lost
 * too much of the answer to cancellation for a double to carry it (the MPFR
 * path answers there). */
bool sinintegral_machine_complex(double are, double aim, double* ore, double* oim);

#endif /* SININTEGRAL_H */
