#ifndef ERFI_H
#define ERFI_H

#include "expr.h"

/* Erfi[z] -- the imaginary error function erfi(z) = erf(i z)/i. */
Expr* builtin_erfi(Expr* res);
void  erfi_init(void);

#endif /* ERFI_H */
