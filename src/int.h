/*
 * int.h
 *
 * Integer-digit / digit-list utilities. Currently:
 *   - IntegerDigits[n], IntegerDigits[n, b], IntegerDigits[n, b, len]
 *   - IntegerLength[n], IntegerLength[n, b]
 *   - DigitCount[n], DigitCount[n, b], DigitCount[n, b, d]
 */

#ifndef MATHILDA_INT_H
#define MATHILDA_INT_H

#include "expr.h"

Expr* builtin_integerdigits(Expr* res);
Expr* builtin_integerlength(Expr* res);
Expr* builtin_digitcount(Expr* res);

void int_init(void);

#endif /* MATHILDA_INT_H */
