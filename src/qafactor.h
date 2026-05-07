/* qafactor.h — Trager algebraic-factoring helpers (Phase G).
 *
 * This module sits on top of qa.{c,h} (Q(α) elements) and qaupoly.{c,h}
 * (Q(α)[x] polynomials) and connects them to picocas's existing
 * polynomial machinery (Resultant, Factor, ZUPoly).
 *
 * Phase G3 — Norm via resultant.
 *
 *   Given f(x, α) ∈ Q(α)[x] and minimal polynomial P_α(y) ∈ Q[y], the
 *   field-norm
 *
 *      N(f)(x) = Resultant_y( P_α(y), g(x, y) )
 *
 *   where g(x, y) ∈ Q[x, y] is f with α textually replaced by y, lies
 *   in Q[x].  When f is irreducible over Q(α), N(f) is a power of an
 *   irreducible Q-factor up to (separable) shifts — this is the
 *   foundation of Trager's algorithm.
 *
 * Implementation strategy.  Rather than reinvent multivariate
 * resultants we serialise the QAUPoly to a picocas Expr (a polynomial
 * in two free symbols, x_name and y_name) and call the existing
 * `internal_resultant`, which dispatches to the Sylvester-matrix
 * routine in poly.c.  The result comes back as a univariate Expr in
 * x_name, expanded into canonical form. */

#ifndef PICOCAS_QAFACTOR_H
#define PICOCAS_QAFACTOR_H

#include "qa.h"
#include "qaupoly.h"

struct Expr;

/* Build P_α(y) as a picocas Expr (a polynomial in `y_name`).
 * Caller owns the returned Expr. */
struct Expr* qaext_to_expr(const QAExt* ext, const char* y_name);

/* Build f(x, α) → g(x, y) ∈ Q[x, y] as a picocas Expr by textual
 * substitution α → y.  The result is a polynomial in two free
 * symbols `x_name` and `y_name`.  Caller owns the returned Expr. */
struct Expr* qaupoly_to_expr(const QAUPoly* f,
                             const char* x_name,
                             const char* y_name);

/* Compute the field norm
 *
 *     N(f)(x) = Resultant_y(P_α(y), g(x, y))   ∈ Q[x]
 *
 * via picocas's own Resultant builtin.  The result is returned as a
 * picocas Expr (a univariate polynomial in `x_name`, post-Expand).
 * Returns NULL if f is the zero polynomial.  Caller owns the returned
 * Expr. */
struct Expr* qaupoly_norm(const QAUPoly* f,
                          const char* x_name,
                          const char* y_name);

#endif
