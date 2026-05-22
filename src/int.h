/*
 * int.h
 *
 * Integer-digit / digit-list utilities. Currently:
 *   - IntegerDigits[n], IntegerDigits[n, b], IntegerDigits[n, b, len]
 */

#ifndef MATHILDA_INT_H
#define MATHILDA_INT_H

#include "expr.h"

Expr* builtin_integerdigits(Expr* res);

void int_init(void);

#endif /* MATHILDA_INT_H */
