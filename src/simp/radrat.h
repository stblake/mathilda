#ifndef MATHILDA_RADRAT_H
#define MATHILDA_RADRAT_H

#include "expr.h"
#include "simp.h"

/*
 * radrat.c -- Fast algebraic normal form for rational functions of
 * MULTIPLE radical generators.
 *
 * Mathilda's Together/Cancel already combine a single radical generator
 * B^(1/q) with its mixed integer/fractional powers (rat.c + the
 * qa_cancel_with_poly_radical path in poly/qafactor.c).  They give up,
 * however, when two or more DISTINCT radical bases occur and the
 * simplification needs the relations between them -- e.g.
 *
 *     D[ Integrate[1/(x^3 (a+b x)^(1/3)), x], x ]
 *
 * is a sum of fractions over a^(k/3) and (a+b x)^(k/3) that only collapse
 * back to the integrand 1/(x^3 (a+b x)^(1/3)) after reducing modulo the
 * generator relations s^3 = a, t^3 = a+b x.
 *
 * simp_radical_rational works in the quotient ring
 *
 *     K[ g_1, ..., g_n, free_vars ] / < g_k^{q_k} - base_k >
 *
 * for the distinct radical bases {base_k} (each provably positive, so the
 * substitution base_k^(p/q) -> g_k^p is branch-safe).  It substitutes the
 * bases to fresh polynomial generators, combines over a common denominator
 * (Together), reduces numerator and denominator modulo the relation ideal
 * (PolynomialRemainder), rationalises the denominator (PolynomialExtended-
 * GCD), substitutes the radicals back, and finishes with a single-radical
 * Cancel.  The result is returned only when its score is STRICTLY smaller
 * than the input's, so the pass can never regress a case.
 *
 * Contract (mirrors simp_trig_rational):
 *   - Input is borrowed (not freed).
 *   - On success: returns a freshly-allocated Expr* the caller must free.
 *   - On no-op / non-applicable / no-improvement: returns NULL, and the
 *     caller (simp_dispatch) continues with the normal pipeline unchanged.
 *
 * Scope: invoked only from Simplify's dispatcher; Together/Cancel keep
 * their existing single-generator behaviour.
 */
Expr* simp_radical_rational(const Expr* input,
                            const AssumeCtx* ctx,
                            const Expr* complexity_func);

#endif /* MATHILDA_RADRAT_H */
