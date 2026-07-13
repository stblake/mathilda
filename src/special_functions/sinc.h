#ifndef SINC_H
#define SINC_H

#include "expr.h"

/* Sinc[z] -- the (unnormalized) cardinal sine  Sin[z]/z, with Sinc[0] = 1.
 * An entire, even function; the removable singularity at 0 is filled in. */
Expr* builtin_sinc(Expr* res);
void  sinc_init(void);

#endif /* SINC_H */
