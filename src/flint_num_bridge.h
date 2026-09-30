/*
 * flint_num_bridge.h
 * ------------------
 * Boundary between Mathilda's numeric Expr scalars and FLINT's arbitrary-
 * precision complex ball arithmetic (arb/acb), for the transcendental special
 * functions. FLINT's acb_hypgeom / acb_dirichlet kernels are rigorous (correct
 * rounding, full complex plane, standard branch cuts), so they back the FLINT`
 * context wrappers exposed at the REPL and, over time, transparent numeric
 * fast paths for the hand-rolled MPFR kernels.
 *
 * Each entry point takes numeric Expr arguments (Integer / BigInt / Rational /
 * Real / MPFR / Complex of those) and returns a freshly-owned numeric Expr
 * (EXPR_MPFR, or Complex[...] when the imaginary part is non-zero), using the
 * working precision implied by the arguments (machine 53-bit floor). It returns
 * NULL — leaving the caller's FLINT`f[...] unevaluated — when an argument is
 * non-numeric / symbolic, or the result is a pole / non-finite, or FLINT is
 * absent (USE_FLINT undefined). No argument is mutated.
 */
#ifndef FLINT_NUM_BRIDGE_H
#define FLINT_NUM_BRIDGE_H

#include "expr.h"

#ifdef __cplusplus
extern "C" {
#endif

Expr* flint_num_zeta(const Expr* s);                       /* Zeta[s]            */
Expr* flint_num_hurwitz_zeta(const Expr* s, const Expr* a);/* HurwitzZeta[s, a]  */
Expr* flint_num_polygamma(const Expr* n, const Expr* z);   /* PolyGamma[n, z]    */
Expr* flint_num_stieltjes(const Expr* n, const Expr* a);   /* StieltjesGamma[n,a]*/

/* The Legendre elliptic integrals, in the PARAMETER convention m = k^2 that
 * Mathematica uses (not the modulus k) -- which is also Arb's, so no argument
 * translation happens here:
 *
 *   EllipticK[m]           = F(Pi/2 | m)
 *   EllipticE[m]           = E(Pi/2 | m)
 *   EllipticF[phi, m]      = Int_0^phi dt / Sqrt(1 - m Sin[t]^2)
 *   EllipticE[phi, m]      = Int_0^phi Sqrt(1 - m Sin[t]^2) dt
 *   EllipticPi[n, m]       = Pi(n; Pi/2 | m)
 *   EllipticPi[n, phi, m]  = Int_0^phi dt / ((1 - n Sin[t]^2) Sqrt(1 - m Sin[t]^2))
 *
 * Arb defines the incomplete forms on -Pi/2 <= Re phi <= Pi/2 and extends them
 * by the quasi-period (F(phi + k Pi | m) = F(phi | m) + 2 k K(m), and likewise
 * for E and Pi), so a phi of any size or complexity is handled here rather than
 * by the caller. EllipticPi with n > 1 crosses the pole at Sin[t]^2 = 1/n and is
 * the Cauchy principal value, which acb_elliptic_pi_inc computes -- a complex
 * value on the principal branch, matching Mathematica. */
Expr* flint_num_elliptic_k(const Expr* m);                              /* EllipticK[m]        */
Expr* flint_num_elliptic_e(const Expr* m);                              /* EllipticE[m]        */
Expr* flint_num_elliptic_f(const Expr* phi, const Expr* m);             /* EllipticF[phi, m]   */
Expr* flint_num_elliptic_e_inc(const Expr* phi, const Expr* m);         /* EllipticE[phi, m]   */
Expr* flint_num_elliptic_pi(const Expr* n, const Expr* m);              /* EllipticPi[n, m]    */
Expr* flint_num_elliptic_pi_inc(const Expr* n, const Expr* phi,
                                const Expr* m);                          /* EllipticPi[n,phi,m] */

/* Registers the FLINT` context numeric builtins (Zeta / HurwitzZeta /
 * PolyGamma / StieltjesGamma). Called from core_init(). No-op without FLINT. */
void flint_num_bridge_init(void);

#ifdef __cplusplus
}
#endif

#endif /* FLINT_NUM_BRIDGE_H */
