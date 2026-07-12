/* risch_field.h — differential-field primitives for the Bronstein Risch algorithm.
 *
 * A monomial extension k(t) of a differential field (k, D) carries a derivation
 * D that acts on the integration variable and on a finite set of tower monomial
 * variables t_1, ..., t_n, each with a known derivative Dt_i (Bronstein,
 * *Symbolic Integration I*, 2nd ed., Ch. 3).  This module provides the shared
 * primitives every downstream Risch phase (canonical representation, Hermite
 * reduction, the residue criterion, the Risch DE) is built on:
 *
 *   - the monomial derivation D[p] = sum_i (Dt_i) d p / d t_i ,
 *   - gcd and exact division in k[t] treating the coefficient ring k = C(other
 *     variables) as a FIELD (so pure-x factors are units — the operations are
 *     over Q(x)[t], not Q[x, t]),
 *   - the normal / special classification of a polynomial (Def. 3.4.2).
 *
 * Arithmetic is grounded in Mathilda's existing Expr builtins (D, Together,
 * PolynomialGCD, Cancel, Coefficient, CoefficientList), which carry their own
 * FLINT fast paths for rational operands.  There is no separate polynomial-list
 * representation: a k[t] element is an Expr that is a polynomial in the monomial
 * variable t with coefficients rational in the remaining variables.
 *
 * Memory: every returned Expr* is freshly allocated and owned by the caller.
 */

#ifndef MATHILDA_RISCH_FIELD_H
#define MATHILDA_RISCH_FIELD_H

#include <stdbool.h>
#include <stddef.h>
#include "expr.h"

/* A monomial derivation, given as a set of (variable, derivative) pairs.  D acts
 * as D[p] = sum_i dvars[i] * (d p / d vars[i]).  The list MUST include the base
 * integration variable (with derivative 1) and every monomial variable with its
 * derivative; any symbol not listed is treated as a constant (derivative 0).
 * The `vars`/`dvars` arrays are owned by this struct; the Expr* elements they
 * point at are BORROWED (typically subexpressions of a caller-held rule list). */
typedef struct {
    size_t nvars;
    Expr** vars;   /* borrowed element pointers */
    Expr** dvars;  /* borrowed element pointers */
} RischDeriv;

/* Populate `out` from a List[Rule[var, Dvar], ...] Expr.  Borrows the rule
 * subexpressions (valid while `rules` lives).  Returns false on malformed input
 * (non-List head, a non-Rule element, or a Rule of arity != 2). */
bool risch_deriv_from_rules(const Expr* rules, RischDeriv* out);

/* Release the vars/dvars arrays allocated by risch_deriv_from_rules.  Does NOT
 * free the borrowed Expr elements. */
void risch_deriv_free(RischDeriv* d);

/* The monomial derivation D[p].  Owned result (evaluated). */
Expr* risch_field_deriv(const Expr* p, const RischDeriv* d);

/* Degree of p in the monomial variable t, or -1 if p is (structurally) zero.
 * Requires p to be a polynomial in t. */
long risch_field_degree_t(const Expr* p, const Expr* t);

/* Monic-in-t gcd of a and b over the field k = C(other variables).  Returns the
 * integer 1 when a and b share only a unit (gcd of degree 0 in t).  Owned. */
Expr* risch_field_gcd_t(const Expr* a, const Expr* b, const Expr* t);

/* Exact quotient a/b in k(t), requiring b | a with a polynomial-in-t result.
 * Owned result, or NULL when the quotient is not a polynomial in t. */
Expr* risch_field_divexact_t(const Expr* a, const Expr* b, const Expr* t);

/* Numerator a and denominator d of f as a rational function in t over k, with d
 * made monic in t (a/d = f).  Both owned. */
void risch_field_num_den_t(const Expr* f, const Expr* t, Expr** a, Expr** d);

/* Polynomial division in k[t]: a = q b + r with deg_t(r) < deg_t(b).  Returns
 * false (and leaves *q, *r untouched) if b is the zero polynomial.  Owned. */
bool risch_field_divmod_t(const Expr* a, const Expr* b, const Expr* t,
                          Expr** q, Expr** r);

/* Extended gcd in k[t]: returns g, u, v (all owned) with u a + v b = g. */
void risch_field_xgcd_t(const Expr* a, const Expr* b, const Expr* t,
                        Expr** g, Expr** u, Expr** v);

/* Diophantine solve for coprime dn, ds: find b, c with b dn + c ds = r and
 * deg_t(b) < deg_t(ds).  Requires gcd_t(dn, ds) = 1.  b, c owned. */
void risch_field_diophantine_t(const Expr* dn, const Expr* ds, const Expr* r,
                               const Expr* t, Expr** b, Expr** c);

/* Def. 3.4.2: p is normal iff gcd(p, D[p]) = 1; special iff p | D[p]. */
bool risch_field_is_normal(const Expr* p, const Expr* t, const RischDeriv* d);
bool risch_field_is_special(const Expr* p, const Expr* t, const RischDeriv* d);

#endif /* MATHILDA_RISCH_FIELD_H */
