/*
 * real.h
 *
 * RealDigits builtin: positional-notation digit expansion of approximate
 * and exact real numbers.
 *
 *   RealDigits[x]                {digits, exp} in base 10, len from precision
 *   RealDigits[x, b]             same but in base b
 *   RealDigits[x, b, len]        len digits
 *   RealDigits[x, b, len, n]     len digits, first one coefficient of b^n
 *
 * The result is { digits-list, exp }.  The first element of `digits-list`
 * is the coefficient of b^(exp-1).  Sign of x is discarded.  For exact
 * rationals with non-terminating expansions, the digit list ends in a
 * nested list of the recurring block.  For inexact reals, digits past the
 * available precision are filled in as the symbol Indeterminate.
 */

#ifndef MATHILDA_REAL_H
#define MATHILDA_REAL_H

#include "expr.h"

Expr* builtin_realdigits(Expr* res);

void real_init(void);

#endif /* MATHILDA_REAL_H */
