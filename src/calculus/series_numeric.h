/*
 * series_numeric.h -- the Taylor series of an INEXACT analytic expression,
 * computed numerically instead of through the exact round trip.
 *
 * Series of an input that carries machine or MPFR numbers is otherwise
 * computed by rationalising the input, expanding it exactly and numericalising
 * the result (internal_rationalize_then_numericalize).  That is right for the
 * shapes the exact engine needs exact coefficients for (leading-term
 * cancellation in a Laurent or Puiseux expansion, PolynomialQuotient, Apart),
 * but for an ANALYTIC expansion it pays for symbolic arithmetic the answer then
 * discards: 2.^(z/2 + 1) rationalises to 2^(z/2 + 1), whose exact series
 * carries Log[2]^k symbolically, and the 30 powers of that series a Schroeder
 * matrix needs cost 0.8 s where the numeric convolutions cost a millisecond.
 *
 * series_inexact_taylor computes the same expansion -- of the same rationalised
 * input, so 0.1 still means 1/10 -- in MPFR at the answer's precision plus
 * guard bits, with a running magnitude bound per coefficient (the "majorant":
 * the size of the largest term that went into it).  A coefficient below its
 * bound's rounding level is an exact zero, as the exact expansion would have
 * found it; one that is neither clearly zero nor resolved to the answer's
 * precision, or any operation the engine does not cover (a pole, a branch
 * point, a symbolic parameter, a complex constant), declines to the exact
 * round trip.  The answer is then the exact path's to the last bit except at a
 * rounding tie closer than 2^-64.
 */
#ifndef SERIES_NUMERIC_H
#define SERIES_NUMERIC_H

#include "expr.h"

/* `res` is a Series[f, {x, x0, n}] call whose arguments contain an inexact
 * number.  Returns the SeriesData the exact round trip would return, or NULL
 * (touching nothing) when the expansion is not a plain Taylor expansion of
 * real analytic operations the engine covers. */
Expr* series_inexact_taylor(const Expr* res);

#endif
