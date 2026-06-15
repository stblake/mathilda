#ifndef AIRYBI_H
#define AIRYBI_H

#include "expr.h"

/* AiryBi[z] -- the Airy function Bi(z), the solution of y'' = z y that grows
 * exponentially as z -> +Infinity. Bi is an entire function of z.
 * AiryBiPrime[z] -- its derivative (minimal head: exact value at 0 only). */
Expr* builtin_airybi(Expr* res);
Expr* builtin_airybiprime(Expr* res);
void  airybi_init(void);

#endif /* AIRYBI_H */
