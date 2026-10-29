/*
 * flint_polymat.h
 * ---------------
 * Polynomial matrices over K[x] for a coefficient field K in {Q, GF(p),
 * Q(theta)}, through FLINT's generic rings (gr / gr_poly / gr_mat), exposed
 * in the FLINT` context:
 *
 *   FLINT`PolynomialHermiteForm[M, x, opts]       Hermite normal form
 *   FLINT`PolynomialDetAdjugate[M, x, opts]       {det, adjugate}
 *   FLINT`PolynomialWeakPopovForm[M, x, s, opts]  shifted weak Popov form
 *   FLINT`PolynomialQuotientRemainder[a, b, x, opts]   {q, r}
 *   FLINT`PolynomialExtendedGCD[a, b, x, opts]    {g, {s, t}}
 *   FLINT`PolynomialFieldGCD[a, b, x, opts]       monic gcd (modular over Q(theta))
 *
 * opts: Modulus -> p (GF(p), p a word-size prime) or Extension -> {t, mu}
 * (Q(theta) = Q[t]/(mu), mu irreducible over Q, elements written as
 * polynomials in t), and for the Hermite form "DeterminantMultiple" -> D.
 * Every routine leaves the call unevaluated when its input is out of scope.
 * The Trager integration package (src/internal/trager/) is the main client:
 * these are the kernels of its divisor-class arithmetic.
 */
#ifndef FLINT_POLYMAT_H
#define FLINT_POLYMAT_H

void flint_polymat_init(void);

#endif /* FLINT_POLYMAT_H */
