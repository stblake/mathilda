/* Mathilda -- AiryAiZero, AiryAiPrimeZero, AiryBiZero, AiryBiPrimeZero: the
 * zeros of the Airy functions. See airyzero.c. */
#ifndef MATHILDA_AIRYZERO_H
#define MATHILDA_AIRYZERO_H

#include "expr.h"

Expr* builtin_airyzero(Expr* res);
/* Register the four heads, their attributes and N route. Called from
 * core_init(). */
void airyzero_init(void);

#endif /* MATHILDA_AIRYZERO_H */
