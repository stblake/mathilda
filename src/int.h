/*
 * int.h
 *
 * Integer-digit / digit-list utilities. Currently:
 *   - IntegerDigits[n], IntegerDigits[n, b], IntegerDigits[n, b, len]
 *   - IntegerLength[n], IntegerLength[n, b]
 *   - IntegerExponent[n], IntegerExponent[n, b]
 *   - DigitCount[n], DigitCount[n, b], DigitCount[n, b, d]
 *   - FromDigits[list], FromDigits[list, b],
 *     FromDigits["string"], FromDigits["string", b]
 *   - IntegerString[n], IntegerString[n, b], IntegerString[n, b, len]
 */

#ifndef MATHILDA_INT_H
#define MATHILDA_INT_H

#include "expr.h"

Expr* builtin_integerdigits(Expr* res);
Expr* builtin_integerlength(Expr* res);
Expr* builtin_integerexponent(Expr* res);
Expr* builtin_digitcount(Expr* res);
Expr* builtin_fromdigits(Expr* res);
Expr* builtin_integerstring(Expr* res);

void int_init(void);

#endif /* MATHILDA_INT_H */
