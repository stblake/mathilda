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

/*
 * MantissaExponent[x]       -> {m, e} such that x = m * 10^e, 1/10 <= |m| < 1.
 * MantissaExponent[x, b]    -> base-b mantissa/exponent.
 *
 * Accepts Integer, BigInt, Rational, Real (machine), and (with USE_MPFR)
 * MPFR-precision real arguments.  Complex arguments emit `::realx`; non-
 * integer bases leave the call unevaluated (we currently only handle
 * integer bases >= 2).
 */
Expr* builtin_mantissa_exponent(Expr* res);

/*
 * RealExponent[x]           -> Log[10, |x|]   (as an inexact number).
 * RealExponent[x, b]        -> Log[b,  |x|].
 *
 * Accepts Integer, BigInt, Rational, Real (machine), MPFR, and any
 * symbolic argument that numericalizes to a real value (Pi, E, ...).
 * The base may be any positive numeric > 1 (Integer, Rational, Real,
 * MPFR, or symbolic like E).  Returns a machine Real unless one of the
 * inputs carries MPFR precision, in which case the result is MPFR at
 * the higher of the input precisions.
 *
 * For zero inputs we follow Mathilda's Accuracy convention -- every
 * zero has Accuracy = Infinity -- so RealExponent[0] = -Infinity.
 *
 * Diagnostics:
 *   - argc not in [1, 2] -> ::argt
 *   - Complex x          -> ::realx
 *   - base <= 1          -> ::ibase
 */
Expr* builtin_real_exponent(Expr* res);

void real_init(void);

#endif /* MATHILDA_REAL_H */
