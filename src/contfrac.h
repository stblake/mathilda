/*
 * contfrac.h — ContinuedFraction builtin.
 *
 * Registers ContinuedFraction[x] / ContinuedFraction[x, n], computing the
 * simple continued-fraction expansion of x. See contfrac.c for the
 * algorithm details.
 */
#ifndef CONTFRAC_H
#define CONTFRAC_H

#include "expr.h"

/* Register the ContinuedFraction / FromContinuedFraction builtins,
 * attributes and docstring hooks. */
void contfrac_init(void);

/* Builtin entry points (exposed for testing / direct dispatch). */
Expr* builtin_continued_fraction(Expr* res);
Expr* builtin_from_continued_fraction(Expr* res);

#endif /* CONTFRAC_H */
