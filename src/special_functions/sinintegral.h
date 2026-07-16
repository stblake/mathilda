#ifndef SININTEGRAL_H
#define SININTEGRAL_H

#include "expr.h"

/* SinIntegral[z] -- the sine integral  Si(z) = Int_0^z Sin[t]/t dt.
 * An entire, odd function with no branch cuts. */
Expr* builtin_sinintegral(Expr* res);
void  sinintegral_init(void);

#endif /* SININTEGRAL_H */
