/*
 * ff_xgcd.h -- the extended gcd over K[u], K = Q[t]/(mu) or Q, by the certified
 * multi-modular driver (ff_modular.h); ffield library, Phase 4.5.
 *
 * WHY.  The Euclidean remainder sequence over Q or a number field swells: two
 * coprime polynomials of degree 26 over Q(omega) took 23 s in FLINT's
 * gr_poly_xgcd (G-perf-2, Trager's kExtGCD), where the cofactors themselves
 * are small.  Images over F_p[t]/(f) cost microseconds.
 *
 * The answer.  g monic, s a + t b = g, deg s < deg b - deg g (s = 0 when
 * deg b = deg g) and deg t < deg a - deg g.  Under these bounds (s, t) is
 * unique: if (s', t') is another pair, (s - s')(a/g) = (t' - t)(b/g) with a/g,
 * b/g coprime, so b/g divides s - s', of lower degree, and s = s'.  So the
 * candidate is the answer of any exact extended gcd normalised the same way
 * (FLINT's).
 *
 * The certificate (exact, in K[u]): g monic of the profile's degree, g | a,
 * g | b, s a + t b = g, and the degree bounds.  Then g divides every common
 * divisor (it is a combination of a and b) and divides both, so g = gcd(a, b).
 *
 * Declined (returns 0): a or b zero; deg a = deg b = deg g (a and b constant
 * multiples of g: the cofactors are not unique); a prime budget exhausted.
 */
#ifndef FF_XGCD_H
#define FF_XGCD_H

#ifdef USE_FLINT
#include "ff_modular.h"

/* g, s, t initialised by the caller.  1 with the certified answer, 0 to
 * decline.  mu monic irreducible over Q, or NULL for Q. */
int ff_xgcd_modular(ff_qtpoly_struct* g, ff_qtpoly_struct* s, ff_qtpoly_struct* t,
                    const ff_qtpoly_struct* a, const ff_qtpoly_struct* b,
                    const fmpq_poly_struct* mu, ff_mm_stats* stats);

#endif /* USE_FLINT */
#endif /* FF_XGCD_H */
