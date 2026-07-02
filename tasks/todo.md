# Simplify numericalisation bug + FLINT speedup

## Problem (two asks)
1. `{0.982981, r} // Simplify` (and minimal `{1., 1} // Simplify -> {1.0, 1.0}`)
   numericalised exact numbers that merely shared a list with an inexact one.
2. Simplify of the exact Goursat antiderivative took ~16 s — could FLINT help?

## Root causes
1. Simplify is not Listable, so a List arg went whole through
   `internal_rationalize_then_numericalize`, whose final step blanket-
   `numericalize`s the entire result tree (rationalize.c:844).
2. `simp_algebraic` (Sqrt->g_i) -> Together/Cancel -> `cancel_recursive`'s
   `PolynomialGCD` over Q[k,x,g1,g2] fell back to the interpreted subresultant
   Euclid (`flint_extension_gcd` rightly declines: no algebraic coeff ring).

## Fixes
- [x] `simp_builtins.c`: thread Simplify over a List when it holds an inexact
      leaf (`simplify_thread_list`), per-element inexact decision, opts preserved.
- [x] `flint_bridge.{c,h}`: `flint_multivariate_gcd_normalized` (monic->WL
      primitive-integer via Gauss's lemma), refactored shared core.
- [x] `rat.c`: `cancel_recursive` tries the normalized FLINT GCD first, falls
      back to the builtin. Scoped to internal GCD — user PolynomialGCD unchanged.

## Review
- `{1., 1} // Simplify == {1., 1}`; nested lists per-element; within-expression
  contagion intact (`Simplify[(x-1.5)/(x^2-2.25)]==1/(1.5+x)`).
- Goursat Simplify: ~16 s -> ~3.7 s (≈4.3×); `FreeQ[Simplify[r], Real] == True`.
- Result unchanged up to a `Log[-1]=iπ` additive constant (valid antiderivative).
- Zero new test failures. Verified green: poly, rat, simp, integrate_goursat,
  zero_test, intrat, apart/together, series, sum, extension, cyclotomic, tower.
- Pre-existing failures (confirmed identical on stashed baseline, NOT ours):
  factor_baseline / factor_phase0 / factor_recombine, parfrac Apart string-match
  (already emitted the fully-factored form), solve hidden-zero radical denest.
