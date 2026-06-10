#ifndef ERFC_H
#define ERFC_H

#include "expr.h"

/* Erfc[z] -- the complementary error function erfc(z) = 1 - erf(z). */
Expr* builtin_erfc(Expr* res);
void  erfc_init(void);

#endif /* ERFC_H */
