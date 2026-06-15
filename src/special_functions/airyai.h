#ifndef AIRYAI_H
#define AIRYAI_H

#include "expr.h"

/* AiryAi[z] -- the Airy function Ai(z), the decaying solution of y'' = z y.
 * AiryAiPrime[z] -- its derivative (minimal head: exact value at 0 only). */
Expr* builtin_airyai(Expr* res);
Expr* builtin_airyaiprime(Expr* res);
void  airyai_init(void);

#endif /* AIRYAI_H */
