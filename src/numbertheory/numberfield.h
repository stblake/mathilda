/*
 * numberfield.{c,h}
 * -----------------
 * Minimal algebraic-number-field layer for the integer Diophantine
 * engines (the Thue solver, and later the elliptic-curve solver): the
 * shared "Tier-3 prerequisite" named in SOLVE_INTEGERS.md.
 *
 * A NumberField models K = Q(theta) with theta a root of a MONIC
 * irreducible integer polynomial f (the defining polynomial).  It carries
 * exactly the data the callers need and can CERTIFY:
 *
 *   - degree n, signature (r1 real embeddings, r2 complex-conjugate
 *     pairs), unit rank r = r1 + r2 - 1;
 *   - the field discriminant d_K;
 *   - a *certified* verdict that the equation order Z[theta] equals the
 *     full ring of integers O_K (the monogenic case) -- this is "Gate 1"
 *     of the Thue solver, decided exactly by Dedekind's criterion.
 *
 * Rigour contract: nf_field_create() returns a field ONLY when Z[theta]
 * is provably maximal (monogenic).  When maximality cannot be certified
 * -- a non-maximal order, a discriminant that will not fully factor, a
 * non-monic / reducible input, or FLINT unavailable -- it returns NULL so
 * the caller DECLINES.  Everything decided here is exact integer / finite-
 * field arithmetic; there is no floating-point in the maximality verdict.
 *
 * Requires FLINT (USE_FLINT).  Without it every constructor returns NULL.
 */
#ifndef NUMBERFIELD_H
#define NUMBERFIELD_H

#include <gmp.h>
#include <stdbool.h>
#include <stdint.h>

typedef struct NumberField NumberField;   /* opaque */

/* Build K from a MONIC integer defining polynomial f given by its
 * coefficients coeffs[0..deg] with coeffs[deg] == 1 and deg >= 2.  Runs
 * Gate 1 (disc, factor, Dedekind at every p with p^2 | disc).  Returns a
 * NumberField on success (Z[theta] certified maximal), else NULL (caller
 * DECLINEs).  Caller owns the result; free with nf_field_free().         */
NumberField* nf_field_create(const mpz_t* coeffs, int deg);

void nf_field_free(NumberField* K);

/* Accessors. */
int  nf_degree(const NumberField* K);
int  nf_r1(const NumberField* K);
int  nf_r2(const NumberField* K);
int  nf_unit_rank(const NumberField* K);
void nf_disc(const NumberField* K, mpz_t out);   /* d_K == disc(f) when monogenic */

/* ------------------------------------------------------------------ *
 *  Gate-1 primitives (exposed for unit testing).                      *
 * ------------------------------------------------------------------ */

/* Signature of the monic irreducible integer polynomial f[0..deg]:
 * sets *r1 to the number of real roots (via an exact Sturm sequence).
 * Returns true on success. */
bool nf_signature(const mpz_t* coeffs, int deg, int* r1);

/* Is the equation order Z[theta] p-maximal for the monic integer poly
 * f[0..deg]?  Decided exactly by Dedekind's criterion.  Returns the
 * verdict; sets *ok=false on an internal failure (e.g. p too large for a
 * single-word modulus, or FLINT absent) in which case the verdict is
 * meaningless and the caller must DECLINE. */
bool nf_dedekind_p_maximal(const mpz_t* coeffs, int deg, const mpz_t p, bool* ok);

/* Whole-field maximality certification: true iff Z[theta] == O_K, decided
 * by factoring disc(f) and running Dedekind at each p with p^2 | disc.
 * Sets *certified=false when the answer cannot be proven (disc will not
 * factor, a prime too large, FLINT absent) -- caller then DECLINEs. */
bool nf_is_maximal_order(const mpz_t* coeffs, int deg, bool* certified);

/* True when the library was built with FLINT and this layer is usable. */
bool nf_available(void);

#endif /* NUMBERFIELD_H */
