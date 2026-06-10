#ifndef EXPINTEGRALEI_H
#define EXPINTEGRALEI_H

#include "expr.h"

/* ExpIntegralEi[z] -- the exponential integral Ei(z), the principal value of
 * -Int_{-z}^Inf e^-t/t dt, with a branch cut along the negative real axis. */
Expr* builtin_expintegralei(Expr* res);
void  expintegralei_init(void);

#endif /* EXPINTEGRALEI_H */
