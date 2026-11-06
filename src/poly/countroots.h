/* Mathilda -- CountRoots[poly, x] / CountRoots[poly, {x, a, b}]: exact root
 * counts of a univariate polynomial. See countroots.c. */
#ifndef MATHILDA_COUNTROOTS_H
#define MATHILDA_COUNTROOTS_H

#include "expr.h"

Expr* builtin_countroots(Expr* res);
/* Register CountRoots and its attributes. Called from core_init(). */
void countroots_init(void);

#endif /* MATHILDA_COUNTROOTS_H */
