/* gbmod.h
 *
 * Self-contained Buchberger Gröbner-basis engine over the prime field
 * GF(p), independent of the characteristic-0 engine in groebner.{c,h}
 * (which carries GMP rational coefficients).  Here coefficients are native
 * `uint64_t` residues in [0, p) with modular arithmetic, so the field is a
 * true GF(p): division is multiplication by the modular inverse.
 *
 * `p` must be PRIME (GF(p) is a field only then; division would otherwise
 * fail).  The caller is responsible for the primality check.  Products stay
 * below p^2, so p must satisfy p < 2^31 for the arithmetic to be exact in
 * `uint64_t`.
 *
 * Only the two monomial orders the GroebnerBasis builtin exposes are
 * supported: GB_ORDER_LEX and GB_ORDER_GREVLEX (shared enum from
 * groebner.h).  Every GFpPoly* is heap-allocated and released with
 * gfp_poly_free; every basis array with gfp_basis_free.
 */
#ifndef GBMOD_H
#define GBMOD_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include "groebner.h"   /* GBOrder, GBPoly */

typedef struct GFpPoly {
    int       n_vars;
    GBOrder   order;      /* GB_ORDER_LEX or GB_ORDER_GREVLEX */
    uint64_t  p;          /* prime modulus */
    int*      exps;       /* n_terms * n_vars, row-major; terms sorted desc  */
    uint64_t* coefs;      /* n_terms, each in [1, p) after normalize         */
    size_t    n_terms;
    size_t    cap;
} GFpPoly;

/* ---- modular arithmetic ---- */
/* Modular inverse of a in GF(p) (p prime, a not a multiple of p). */
uint64_t gfp_inv(uint64_t a, uint64_t p);
/* Trial-division primality test (adequate for p < 2^31). */
bool     gfp_is_prime(uint64_t p);

/* ---- construction / destruction ---- */
GFpPoly* gfp_poly_new(int n_vars, GBOrder order, uint64_t p);
GFpPoly* gfp_poly_copy(const GFpPoly* a);
void     gfp_poly_free(GFpPoly* a);
/* Append a term (unsorted); caller MUST gfp_poly_normalize afterwards. */
void     gfp_poly_push_term(GFpPoly* a, const int* exps, uint64_t coef);
/* Sort desc by order, merge equal monomials (mod p), drop zero coefs. */
void     gfp_poly_normalize(GFpPoly* a);

/* ---- queries ---- */
bool       gfp_poly_is_zero(const GFpPoly* a);
bool       gfp_poly_is_constant(const GFpPoly* a);
const int* gfp_poly_lm(const GFpPoly* a);      /* leading exps; NULL if zero */
uint64_t   gfp_poly_lc(const GFpPoly* a);      /* leading coef; 0 if zero    */
int        gfp_degree_in(const GFpPoly* a, int var);  /* max exp of `var`    */
bool       gfp_contains_var(const GFpPoly* a, int var);

/* ---- conversion ---- */
/* Reduce a rational GBPoly into GF(p) under `order`.  Returns NULL when any
 * coefficient's denominator is divisible by p (no image in GF(p)). */
GFpPoly* gfp_from_gbpoly(const GBPoly* g, GBOrder order, uint64_t p);
/* Lift a GF(p) poly back to a rational GBPoly with integer coefficients in
 * [0, p) (for rendering a modular basis through gb_to_expr). */
GBPoly*  gbpoly_from_gfp(const GFpPoly* a);

/* ---- Gröbner engine ---- */
/* Reduced Gröbner basis of <F> over GF(p).  All inputs share n_vars/order/p.
 * Returns a heap array of GFpPoly* (sorted ascending by leading monomial,
 * Mathematica convention) and sets *out_n; NULL/0 for the zero ideal. */
GFpPoly** gfp_buchberger(GFpPoly* const* F, size_t n, size_t* out_n);
void      gfp_basis_free(GFpPoly** G, size_t n);

/* ---- solving support ---- */
/* Substitute var -> value; returns a fresh, normalized poly. */
GFpPoly* gfp_poly_subst(const GFpPoly* a, int var, uint64_t value);
/* True iff every term is constant or a pure power of `var` (free of others). */
bool     gfp_univariate_in(const GFpPoly* a, int var);
/* Evaluate a poly that is univariate in `var` at var = r (others ignored). */
uint64_t gfp_eval_univariate(const GFpPoly* a, int var, uint64_t r);
/* Evaluate at a complete assignment `values` (length n_vars). */
uint64_t gfp_eval_full(const GFpPoly* a, const uint64_t* values);

#endif /* GBMOD_H */
