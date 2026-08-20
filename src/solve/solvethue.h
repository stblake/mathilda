/*
 * solvethue.{c,h}
 * ---------------
 * Thue-equation engine: the complete integer solution set of an
 * irreducible binary form equation
 *
 *     sum_{j=0}^{n} form[j] * x^(n-j) * y^j  ==  m,     n >= 3,
 *
 * via the Tzanakis-de Weger method (reduce to unit equations in the
 * number field Q(theta), theta a root of F(t,1); Baker's linear forms in
 * logarithms bound the unit exponents; LLL shrinks the bound; enumerate
 * and verify).  Built on the number-field layer (numberfield.{c,h},
 * nfunits.{c,h}).  Requires FLINT.
 *
 * Contract: PROVABLY COMPLETE or DECLINE.  thue_solve_binary_form()
 * returns the complete certified solution set, or signals decline; it
 * never returns an incomplete set as if complete.
 */
#ifndef SOLVETHUE_H
#define SOLVETHUE_H

#include "expr.h"
#include <gmp.h>

typedef struct { mpz_t x, y; } ThueSol;

/* Complete, certified solution set of the Thue equation defined by
 * form[0..n] and m.  On success returns the number of solutions (>= 0)
 * and sets *out to a malloc'd array (caller frees each x,y and the array
 * via thue_sols_free).  Returns -1 to DECLINE (cannot certify: |m| != 1
 * in the current scope, non-monogenic field, uncertified units, degree
 * out of range, or the rigorous exponent bound is unavailable). */
int thue_solve_binary_form(const mpz_t* form, int n, const mpz_t m, ThueSol** out);

/* Test-only: enumerate unit exponents up to an EXPLICIT bound instead of
 * the rigorous Baker bound.  NOT complete unless `bound` provably covers
 * all solutions -- used only to validate the reconstruction mechanism. */
int thue_solve_binary_form_bounded(const mpz_t* form, int n, const mpz_t m,
                                   int bound, ThueSol** out);

void thue_sols_free(ThueSol* sols, int count);

void solvethue_init(void);

#endif /* SOLVETHUE_H */
