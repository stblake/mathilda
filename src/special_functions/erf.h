#ifndef ERF_H
#define ERF_H

#include "expr.h"

/* Erf[z]      -- the error function erf(z).
 * Erf[z0, z1] -- the generalized error function erf(z1) - erf(z0). */
Expr* builtin_erf(Expr* res);
void  erf_init(void);

#endif /* ERF_H */
