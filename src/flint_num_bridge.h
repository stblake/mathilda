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

/* Registers the FLINT` context numeric builtins (Zeta / HurwitzZeta /
 * PolyGamma / StieltjesGamma). Called from core_init(). No-op without FLINT. */
void flint_num_bridge_init(void);

#ifdef __cplusplus
}
#endif

#endif /* FLINT_NUM_BRIDGE_H */
