/*
 * simp_trigexp_zero.h — exact symbolic zero-test for rational functions of a
 * single exponential kernel.
 *
 * Many trigonometric / exponential identities that the general complexity-
 * guided Simplify search cannot prove (and blows up on) become *rational
 * functions of one variable* once every trig/exp kernel is written over
 * t = E^(I x). The canonical example is a Risch antiderivative diff-back
 * `D[G] - f` (SIMPLIFY_GAPS.md Families 1 & 3): a multiple-angle
 * (or symbolic-parameter, I-laden) combination that is identically 0 but
 * defeats `Simplify` — `Simplify[D[Sec^3-antideriv] - Sec^3]` runs >40 s
 * without terminating; the symbolic `1/(a+b Sin x)` diff-back hangs >90 s.
 *
 * The procedure here is exact and terminating, with NO numeric sampling
 * (a decision procedure must never certify by sampling):
 *
 *   1. TrigToExp[e]                       → exponentials of k·(I·x)
 *   2. kernelize  E^(k I x) → t^k         (all integer k, Laurent), in C —
 *      the surface `/. E^(m_. I x):>t^m` rule mis-fires on Optional-in-
 *      Orderless-Times, so the substitution is done structurally here
 *   3. Together                           → one rational P(t,params)/Q(t,params)
 *   4. numerator ≡ 0 ?  (Expand + literal-zero test)
 *
 * Coefficients may be rational functions of other free symbols (symbolic
 * parameters a, b, …). A numerator that Expands to literal 0 proves the
 * identity; a nonzero numerator with only the kernel t free (pure Gaussian-
 * rational coefficients) rigorously disproves it; anything else declines.
 */
#ifndef MATHILDA_SIMP_TRIGEXP_ZERO_H
#define MATHILDA_SIMP_TRIGEXP_ZERO_H

#include "expr.h"

typedef enum {
    TRIGEXP_ZERO_FALSE   = 0,  /* proven NOT identically zero (rigorous)      */
    TRIGEXP_ZERO_TRUE    = 1,  /* proven identically zero (rigorous)          */
    TRIGEXP_ZERO_UNKNOWN = 2   /* declined: not a single-kernel rational form */
} TrigExpZeroResult;

/* Decide whether `e` is identically zero as a rational function of a single
 * exponential kernel t = E^(I·var), auto-detecting `var` as the unique free
 * symbol occurring inside trig/hyperbolic/exponential arguments. `e` is
 * read-only and is NOT consumed. Declines (UNKNOWN) on: no or more than one
 * kernel variable, a non-integer kernel power, or any residual dependence on
 * `var` outside an exponential (bare `var`, `Log[var]`, `Sqrt[var]`, a tower). */
TrigExpZeroResult trigexp_rational_is_zero(const Expr* e);

/* Simplify seed. Returns a fresh Integer 0 iff `trigexp_rational_is_zero(e)`
 * proves `e` identically zero; otherwise NULL (decline — contributes nothing
 * to the candidate set, never a worse form). `e` is read-only. Cheaply gated
 * so non-vanishing / non-trig inputs bail before any Together. */
Expr* transform_trigexp_vanish(const Expr* e);

#endif /* MATHILDA_SIMP_TRIGEXP_ZERO_H */
