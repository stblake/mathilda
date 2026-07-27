#ifndef EXPINTEGRALEI_H
#define EXPINTEGRALEI_H

#include "expr.h"
#include <stdbool.h>

/* ExpIntegralEi[z] -- the exponential integral Ei(z), the principal value of
 * -Int_{-z}^Inf e^-t/t dt, with a branch cut along the negative real axis. */
Expr* builtin_expintegralei(Expr* res);
void  expintegralei_init(void);


/* ExpIntegralEi for a machine complex argument, in the shared kernel ABI.  Returns false
 * to DECLINE — at the singularity, and wherever the ascending series has lost
 * too much of the answer to cancellation for a double to carry it (the MPFR
 * path answers there). */
bool expintegralei_machine_complex(double are, double aim, double* ore, double* oim);

#endif /* EXPINTEGRALEI_H */
