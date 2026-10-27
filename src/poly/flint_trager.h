/*
 * flint_trager.h
 * --------------
 * Factoring over an algebraic number field K, any number of variables, on
 * FLINT: Trager's norm method with the certified modular field GCD of
 * flint_bridge.c.  This is the engine behind Factor[p, Extension -> ...],
 * Factor[p, GaussianIntegers -> True] and the default factorisation over Q(i)
 * of a polynomial with Gaussian coefficients.
 *
 * Without FLINT (USE_FLINT undefined) the entry point exists and declines.
 */
#ifndef FLINT_TRAGER_H
#define FLINT_TRAGER_H

#include <stddef.h>
#include "expr.h"

/*
 * Factor the polynomial `poly` over K = Q(ext[0], ..., ext[n_ext - 1], and the
 * algebraic constants occurring in poly) -- radicals, Root objects and
 * Gaussian coefficients alike; with n_ext == 0 only poly's own constants
 * generate K (Extension -> Automatic).
 *
 * Returns List[c, List[f_1, e_1], ..., List[f_k, e_k]] with poly ==
 * c f_1^e_1 ... f_k^e_k exactly (certified before returning), each f_i
 * irreducible over K and spelled in the caller's own radicals, primitive over
 * Z (over Z[i] when K = Q(i)) in that basis, with its sign (unit) normalised;
 * c is an element of K.  The factors are in no particular order.
 *
 * Returns NULL to DECLINE: no algebraic constant at all (K = Q: the caller's
 * ordinary factoriser applies), a constant poly, a non-polynomial (a negative
 * or symbolic power, a non-algebraic head, an inexact number, a parametric
 * radical such as Sqrt[x]), a field the common-field construction cannot
 * build, or no Trager shift found.  Never mutates its arguments.
 */
Expr* flint_factor_over_field(const Expr* poly, Expr* const* ext, size_t n_ext);

/* The same, as the evaluated product c f_1^e_1 ... f_k^e_k (NULL to decline). */
Expr* flint_factor_over_field_product(const Expr* poly, Expr* const* ext, size_t n_ext);

/* 1 when `e` contains a non-real Gaussian-rational number (Complex[a, b] with
 * b != 0 and a, b exact), else 0.  Cheap structural walk. */
int flint_trager_has_complex(const Expr* e);

#endif /* FLINT_TRAGER_H */
