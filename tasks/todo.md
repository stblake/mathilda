# Simplify trig stress-test gap closure

Plan: `~/.claude/plans/let-s-test-the-implementation-bright-knuth.md`
Corpus: 101 entries; baseline 84/101 → 0.

## Tier A — structural
- [x] **A2** inverse-trig complementary pairs (ArcSin+ArcCos, ArcTan+ArcCot, ArcSec+ArcCsc → π/2).
      Phase-0e pre-pass in `simp_search.c` + `transform_invtrig_complement` in `simp_tan_add.c`.
      Closes #74 (bonus: bare sums → π/2, `+y` too). Corpus 84 → **85**. Tests green.
- [x] **A1** constant-phase (root-of-unity) angle-expansion fallback in `simp_trigexp_zero.c`
      (`expand_pi_phase`), authoritative over a spurious opaque-FALSE. Closes #13, #97, #98
      AND fixes a latent `PossibleZeroQ` wrong-answer bug (returned False for true identities).
      Corpus **88/101**. All regression suites green; soft-fails proven pre-existing; valgrind clean.

## Tier B — algebraic-number
- [x] **B1** trig at rational multiples of π → qqbar. `to_qqbar` + `is_constant_algebraic`
      Sin/Cos/Tan/Cot/Sec/Csc-at-rational-π cases (`trig_arg_pi_rational`, roots of unity);
      top-level constant-algebraic fold in `builtin_simplify` (trig consts bypass simp_search's
      Phase-0d via the SHAPE_TRIG router). Closes #89, #90, #91. Corpus **91/101**. Bonus:
      `RootReduce[Sin[Pi/7]]` now exact. All qqbar-consumer suites green; valgrind clean.
- [x] **B2** inverse-trig rational-π addition. `simp_invtrig_combo_is_zero` (`simp_builtins.c`):
      build `exp(i·e)` from Euler closed forms, prove `== 1` exactly via qqbar (e ≡ 0 mod 2π),
      numeric screen selects branch k=round(e/2π), keep iff k==0. Closes #92, #94. Corpus **93/101**.

## Review (DONE — v0.327)

**Outcome: 84/101 → 93/101. All 9 genuine gaps closed; latent `PossibleZeroQ`
wrong-answer bug (returned False for true affine-phase identities) fixed.**

Remaining 8 nonzero are all correct: typos #64 (`Tan-Sec+1`), #93 (`ArcTan[4/3]`);
conditional #75/76/77/83/100/101 (need `Assumptions`; match Mathematica).

Files: `simp_tan_add.c` (A2 complement transform), `simp_search.c` (A2 Phase-0e
pre-pass), `simp_trigexp_zero.c` (A1 `expand_pi_phase` + authoritative-over-FALSE),
`flint_qqbar.c` (B1 trig-at-rational-π in `to_qqbar`/`is_constant_algebraic`),
`simp_builtins.c` (B1 top-level constant-algebraic fold + B2 `simp_invtrig_combo_is_zero`).

Verification: full corpus re-run; 22 regression suites green; the 4 soft-fails
proven pre-existing via base A/B (identical set); valgrind clean on all tiers
(no leak traces to new code); `make check-messages`/`check-c99` pass. Tests added:
`test_tez_constant_phase`, `test_invtrig_complementary_pairs`,
`test_invtrig_rational_angle_addition`, `test_simplify_trig_rational_pi_products`.
Docs: `docs/spec/builtins/simplification.md` + changelog `2026-10-05.md`; version 0.327.

Not committed/tagged (awaiting user).

## Hyperbolic corpus (DONE — v0.328)

**Outcome: generated the hyperbolic analogue (93 entries, Osborn's rule; exact-
value/Machin sections omitted — no hyperbolic analogue). 85/93 reduce; the 8
remaining are conditional (match Mathematica without `Assumptions`).**

Found 3 of my own construction errors first (Osborn sign flips) via numeric
verification: #12 `Cosh[iπ/2−x]=−i Sinh` (needed +), #13 `Tanh[iπ/2−x]=−Coth`
(needed +), #58 `Sinh⁶+Cosh⁶=(15Cosh2x+Cosh6x)/16` (the (5+3Cos4x)/8 trig form
has no hyperbolic analogue). Then 4 genuine gaps fixed:
- #13/#16/#90 — hyperbolic heads with affine IMAGINARY phase `k x + i c Pi`:
  extended `expand_pi_phase` (`simp_trigexp_zero.c`) with the hyperbolic addition
  formula; `Cosh[icπ]=Cos[cπ]`, `Sinh[icπ]=i Sin[cπ]` fold to root-of-unity coeffs.
- #29 — `Coth[x] - Csch[x] = Tanh[x/2]` reciprocal-difference half-angle rule
  added to `HalfAngle` (`simp_trig_roundtrip.c`); TrigReduce leaves the
  `(Cosh-1)/Sinh` form split as Coth−Csch.

Conditional (correct, left bare; analogues of trig #75/76/77/100/101): #66/#68
(`Sqrt[x-1]Sqrt[x+1]` vs `Sqrt[x²-1]`), #74 (ArcCosh↔Log), #75/76/77 (inverse-hyp
addition), #92 (`ArcTanh[Tanh[x]]`), #93 (`Sqrt[e^{2x}]`; folds under `Reals`).

Verified: trig corpus still 93/101 (no HalfAngle regression); 19 regression
suites green; valgrind clean. Tests: `test_trigexp_zero.c::{test_tez_hyperbolic_phase,
test_tez_hyperbolic_structural}`. Docs + changelog; version 0.328.

## Cross-cutting (per step)
- [ ] Re-run full corpus (scratchpad/corpus.m); 84 baseline must stay green.
- [ ] Regression: `invtrig_simplify_tests simplify_tests trigrat_tests trigexp_zero_tests fullsimplify_corpus_tests`
      + Risch diff-back consumers (zero_test.c, risch_util.c) for A1.
- [ ] Add corpus regression fixture; valgrind new paths; version bump + changelog + docs/spec.

## Notes
- Corpus typos (Mathilda correct, report back): #64 (denom `Tan-Sec+1`), #93 (`ArcTan[4/3]`).
- Correct refusals (no change): #75, #76, #77, #83, #100, #101.
- Pre-existing soft-fail (NOT mine): `test_simplify_algebraic_u_power_extraction` expected
  string `x^2^(3/2)` mis-parses (right-assoc ^); base binary gives identical output.
