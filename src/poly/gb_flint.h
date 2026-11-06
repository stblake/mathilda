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
#include <gmp.h>
#include "groebner.h"

/* Over Q: F rational GBPolys sharing one order. */
int gb_flint_buchberger(GBPoly* const* F, size_t n, GBPoly*** out, size_t* out_n);

/* Over GF(p), p prime: F rational GBPolys, reduced mod p. */
int gb_flint_buchberger_mod(GBPoly* const* F, size_t n, uint64_t p,
                            GBPoly*** out, size_t* out_n);

/* CoefficientDomain -> InexactNumbers: the reduced Lex basis of a
 * zero-dimensional ideal over Q, monic, with dyadic coefficients justified to
 * `bits` bits (DRL basis exactly, the Lex shape from modular FGLM, the
 * coefficients by MPFR linear algebra at two precisions that must agree).
 * Returns 1 with a fresh basis, 0 to decline (use the exact basis). */
int gb_flint_lex_inexact(GBPoly* const* F, size_t n, long bits,
                         GBPoly*** out, size_t* out_n);

/* The quotient ring of a zero-dimensional ideal over Q for the eigenvalue
 * method: d standard monomials of the DRL basis (*B_out, d rows of n_vars
 * ints) and the multiplication matrices (*M_out, n_vars blocks of d x d mpq,
 * M[(v d + r) d + c] = coefficient of b_r in NF(x_v b_c)).  Returns 1 (d = 0
 * for the unit ideal), 0 to decline.  Caller frees both (mpq_clear each). */
int gb_flint_quotient(GBPoly* const* F, size_t n, int maxd,
                      int** B_out, int* d_out, mpq_t** M_out);

/* Lex over Q: the modular shape-position route first (on = 1, the default)
 * or not (0) -- for A/B tests against the exact routes; the environment
 * variable MATHILDA_GB_SHAPE_OFF does the same.  Returns the number of bases
 * the route has produced so far in this process. */
long gb_flint_shape_route(int on);

/* DegreeReverseLexicographic / DegreeLexicographic over Q: the certified
 * multi-modular route first (on = 1, the default) or not (0); the environment
 * variable MATHILDA_GB_MODULAR_OFF does the same.  Returns the number of bases
 * the route has produced so far. */
long gb_flint_modular_route(int on);

#endif /* GB_FLINT_H */
