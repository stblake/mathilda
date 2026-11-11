/* nroots.h — NRoots[lhs == rhs, var, opts]
 *
 * Numerical root-finding for a univariate polynomial equation.  Returns a
 * disjunction of equations  var==r1 || var==r2 || ...  approximating every
 * root (repeating identical equations for roots of multiplicity > 1, per the
 * Wolfram Language).  Works at machine and arbitrary precision, for real and
 * complex coefficients, with three selectable methods: "Aberth" (default),
 * "CompanionMatrix", and "JenkinsTraub".
 *
 * Memory contract: never frees `res`; returns a fresh Expr* or NULL.
 */
#ifndef MATHILDA_NROOTS_H
#define MATHILDA_NROOTS_H

#include "expr.h"
#include <stdbool.h>
#include <stdint.h>

Expr* builtin_nroots(Expr* res);

/* N[Root[f, k]] for a polynomial f with NON-REAL coefficients (the real case
 * is root_numeric.c's).  `poly` is f in the variable `var` (both borrowed).
 * Returns the k-th root (1-based) at `target_bits` (a machine Complex/Real
 * when `want_machine`) in the reference's numbering for such a polynomial:
 * by real part, and among roots of equal real part the non-real ones first
 * (by imaginary part), then the real one.  Roots are counted with
 * multiplicity.  NULL (no message) when a coefficient is not numeric, the
 * polynomial is constant, k is out of range, the solve fails, or MPFR is
 * compiled out. */
Expr* nroots_root_complex(Expr* poly, Expr* var, int64_t k, long target_bits,
                          bool want_machine);
void  nroots_init(void);

#endif /* MATHILDA_NROOTS_H */
