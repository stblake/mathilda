/* Mathilda -- shared helpers for the Chebyshev polynomial builtins.
 *
 * ChebyshevT (first kind, T_n) and ChebyshevU (second kind, U_n) share the same
 * evaluation substrate: the argument-shape predicates, the explicit-polynomial
 * assembler, the MPFR result builder, the output-precision rule, and the value
 * at x = 0 (Cos[n Pi/2], which is identical for both families). Those pieces are
 * collected here so each family file carries only the parts that genuinely
 * differ -- the recurrence seed, the negative-order reflection, and the
 * non-integer-order closed forms.
 *
 * Everything here is pure (no builtin state) and family-agnostic.
 */
#ifndef MATHILDA_CHEBYSHEV_COMMON_H
#define MATHILDA_CHEBYSHEV_COMMON_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <gmp.h>

#include "expr.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Largest order for which the explicit polynomial is built; beyond this a
 * symbolic/exact-x call stays symbolic. */
#define CHEB_POLY_CAP 2000

/* Largest |order| for which the numeric recurrence runs; beyond this an inexact
 * call stays symbolic rather than looping for an absurd number of steps. */
#define CHEB_NUMERIC_DEG_CAP 2000000

/* ------------------------------------------------------------------ */
/* Argument-shape predicates / extractors                             */
/* ------------------------------------------------------------------ */

/* Recognise an exact integer (machine or GMP) argument that fits a long. */
bool cheb_exact_int(const Expr* e, int64_t* out);

/* True if `e` is exactly the integer `v`. */
bool cheb_is_int(const Expr* e, int64_t v);

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
bool cheb_is_inexact_leaf(const Expr* e);

/* True if `e` carries an inexact numeric leaf (also inside Complex[..]). */
bool cheb_is_inexact(const Expr* e);

/* True iff `e` is a pure numeric quantity (Integer/Real/MPFR/BigInt, or
 * Rational/Complex of such), i.e. a numeric series actually collapsed. */
bool cheb_is_number_leaf(const Expr* e);

/* Machine-double value of a real/integer leaf (declines MPFR/BigInt). */
bool cheb_leaf_double(const Expr* e, double* out);

/* If `x_e` is a machine real, or a Complex of machine real/integer parts, fill
 * (*xr, *xi) with its double value and set *real; return false for anything
 * carrying MPFR/BigInt (which must stay on the arbitrary-precision path). */
bool cheb_machine_xy(const Expr* x_e, double* xr, double* xi, bool* real);

/* ------------------------------------------------------------------ */
/* Explicit integer-coefficient polynomial assembly                   */
/* ------------------------------------------------------------------ */

/* Allocate/free an array of `len` mpz_t, each initialised to 0. */
mpz_t* cheb_mpz_array_new(size_t len);
void   cheb_mpz_array_free(mpz_t* a, size_t len);

/* Build the monomial polynomial Sum_i c[i] x^i (i = 0..deg), skipping zero
 * coefficients, then evaluate once. Returns integer 0 if all vanish. */
Expr* cheb_poly_from_coeffs(mpz_t* c, size_t deg, const Expr* x);

/* ------------------------------------------------------------------ */
/* Numeric result helpers (MPFR)                                       */
/* ------------------------------------------------------------------ */

#ifdef USE_MPFR
/* Build a numeric result Expr from an (re, im) MPFR pair at out_prec bits: a
 * real leaf when `force_real` or im rounds to zero, else Complex[..]. */
Expr* cheb_make_result(const mpfr_t re, const mpfr_t im,
                       mpfr_prec_t out_prec, bool force_real);

/* Output precision: the minimum precision among the inexact arguments. */
mpfr_prec_t cheb_out_prec(const Expr* n_e, const Expr* x_e);

/* If n is an inexact real (or MPFR) whose value is an integer within the cap,
 * return true and set *k. An inexact *integer-valued* order evaluates through
 * the stable recurrence rather than the general closed form, which would leave
 * a spurious tiny imaginary part at integer order. */
bool cheb_inexact_int_order(const Expr* n_e, int64_t* k);

/* T_n(x) (kind 0) or U_n(x) (kind 1) for a non-negative integer degree n of
 * any size and a REAL inexact x, by the closed forms
 *   T: cos(n acos x) on [-1, 1], sgn(x)^n cosh(n acosh |x|) off it;
 *   U: sin((n+1) t)/sin t, (+-1)^n (n+1) at x = +-1, sgn(x)^n sinh/sinh off it,
 * computed with log2(n) guard bits so the product n t keeps out_prec bits:
 * O(1) work where the recurrence takes n steps.  The orders past the
 * recurrence's cap (ChebyshevT[2000001, 0.3], ChebyshevU[10^7, 0.5]) and the
 * bignum orders (ChebyshevT[10^30, 0.5`50], correct to all 50 places) go
 * here (ChebyshevT#2, #4, ChebyshevU#4).  NULL for a complex or non-numeric x. */
Expr* cheb_closed_form_real(int kind, const mpz_t n, const Expr* x_e, mpfr_prec_t out_prec);
#endif /* USE_MPFR */

/* T_n(x) (kind 0) or U_n(x) (kind 1) at an exact rational x for a degree past
 * the explicit-polynomial cap: [T_n, T_(n-1)] = M^(n-1) [x, 1] (U: [2x, 1]) for
 * M = {{2x, -1}, {1, 0}}, by binary powering of the integer matrix
 * q M = {{2p, -q}, {q, 0}} (x = p/q) -- O(log n) products, the value exact
 * (ChebyshevT[2001, 1/2] = -1, ChebyshevT#3, ChebyshevU#3).  NULL for any
 * other x, or when the result would exceed ~2^24 bits. */
Expr* cheb_exact_rational(int kind, unsigned long n, const Expr* x_e);

/* ------------------------------------------------------------------ */
/* Shared special value                                                */
/* ------------------------------------------------------------------ */

/* T_n(0) = U_n(0) = Cos[n Pi/2]. Built as an expression and evaluated, so it
 * serves a symbolic n (-> Cos[n Pi/2]) and an inexact n (-> a number) alike. */
Expr* cheb_value_at_zero(const Expr* n_e);

#endif /* MATHILDA_CHEBYSHEV_COMMON_H */
