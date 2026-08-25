/*
 * reduce_realdiag.h
 *
 * General univariate real sign diagram for `Reduce` (Phase 9): the fallback used
 * when the exact polynomial engine (reduce_univar.c) declines because an atom is
 * a real radical / rational-pole / bounded-domain-transcendental expression
 * rather than a polynomial.
 *
 * Method (the same 1-D cell decomposition, generalised):
 *   1. breakpoints = the real roots (via Solve over Reals, soft) of every atom
 *      numerator and denominator, plus the domain boundaries of every
 *      partial-domain node (radicand==0, Log arg==0, ArcSin arg==±1, ...);
 *   2. sort/dedup with an exact (qqbar) compare that falls back to a numeric
 *      sign when the constant is transcendental (e.g. a multiple of Pi);
 *   3. at one sample per cell, evaluate the formula's truth under a DOMAIN GATE
 *      (a point where any real-domain constraint fails, or a denominator
 *      vanishes, is excluded) — equations/Unequal decided by evaluate() (so a
 *      radical/transcendental identity resolves via PossibleZeroQ), orderings by
 *      an exact-then-numeric sign;
 *   4. emit the union of satisfying cells (shared rru_emit_sign_diagram).
 *
 * Soundness is preserved: a free parameter, an undecidable sign, or an
 * unsupported domain node makes the engine return NULL and Reduce stays
 * unevaluated — never a guessed formula.
 */
#ifndef REDUCE_REALDIAG_H
#define REDUCE_REALDIAG_H

#include "expr.h"
#include "reduce_form.h"
#include "reduce_opts.h"

/* Solve the DNF formula `F` in the single real variable `var` over the reals,
 * tolerating radical / rational-pole / bounded-domain-transcendental atoms.
 * `opts` supplies WorkingPrecision (the numeric-fallback tolerance for
 * transcendental sign decisions) and forwards Cubics / Quartics onto the soft
 * root isolation.  Returns a freshly-owned Expr (True / False / a logical
 * combination), or NULL to decline. */
Expr* reduce_univar_general(const RForm* F, const Expr* var, Expr** vars, int nv,
                            const ReduceOpts* opts);

#endif /* REDUCE_REALDIAG_H */
