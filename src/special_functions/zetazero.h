/* Mathilda -- ZetaZero[k]: the k-th nontrivial zero of Zeta. See zetazero.c. */
#ifndef MATHILDA_ZETAZERO_H
#define MATHILDA_ZETAZERO_H

#include "expr.h"

Expr* builtin_zetazero(Expr* res);
/* Register ZetaZero, its attributes and N route. Called from core_init(). */
void zetazero_init(void);

#endif /* MATHILDA_ZETAZERO_H */
