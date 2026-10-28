/* gb_flint.h -- Buchberger's algorithm on FLINT arithmetic (gb_flint.c).
 *
 * gb_buchberger, gb_groebner_walk and gfp_buchberger try these first.  Both
 * return 1 with the reduced basis in *out / *out_n -- the form those engines
 * return: primitive over Z with a positive leading coefficient (over GF(p):
 * monic, coefficients 0..p-1), sorted ascending by leading monomial, {1} for
 * the unit ideal, empty for the zero ideal -- and 0 to decline: an order with
 * no FLINT ordering (an elimination block order, a weight matrix other than
 * DegreeLexicographic's), a denominator divisible by p, USE_FLINT=0, or
 * MATHILDA_GB_FLINT=off.  The inputs are read, never taken. */
#ifndef GB_FLINT_H
#define GB_FLINT_H

#include <stddef.h>
#include <stdint.h>
#include "groebner.h"

/* Over Q: F rational GBPolys sharing one order. */
int gb_flint_buchberger(GBPoly* const* F, size_t n, GBPoly*** out, size_t* out_n);

/* Over GF(p), p prime: F rational GBPolys, reduced mod p. */
int gb_flint_buchberger_mod(GBPoly* const* F, size_t n, uint64_t p,
                            GBPoly*** out, size_t* out_n);

#endif /* GB_FLINT_H */
