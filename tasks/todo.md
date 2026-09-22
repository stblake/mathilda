# TODO: Fix Mathilda↔Mathematica divergences (MATHILDA_DIVERGENCES.md)

Plan: `/Users/user/.claude/plans/let-s-improve-mathilda-per-zazzy-creek.md`

Fix the 10 A-series correctness bugs + B1/B2/B3, add FLINT nmod_poly fast path for
the GF(p) polynomial ops (section C), then remove the `.m` workarounds and re-verify.
Every fix: extensive unit tests, no leaks, packed/NDArray + Compile surfaces respected.

## Phase 1 — A1 packed/NDArray iterator lists (silent-wrong; top priority)
- [x] Fix `iter_spec_parse` (src/iter.c:98) — materialise EXPR_NDARRAY bound so packed lists iterate (fixes Do/Table/Sum/Product/plotters). VERIFIED.
- [x] `Scan` (src/funcprog.c:1314) — NOT a bug; has correct NDArray fast path (iterates per-element). Earlier "whole buffer" reading was a Reap-grouping misread.
- [ ] Fix `Manipulate`/`Animate` discrete controls (src/graphics/manipulate.c:139)
- [ ] Tests: packed lists through Do/Table/Sum/Product/Scan (test_iter.c / test_packed_list.c)
- [ ] Audits green: check-packed-aware, check-nd-surfaces, check-fastpath-sweep

## Phase 2 — Polynomial correctness (A6→A2, A3, A7) — DONE
- [x] A6: coeff_mod_int/rational_mod_int in polynomial_mod_single reduce Rational + BigInt coeffs (poly.c). VERIFIED.
- [x] A2: PolynomialExtendedGCD modular now correct (FLINT xgcd + classical). VERIFIED {6+x,{6,1}}.
- [x] A3: Modulus for PolynomialQuotient/Remainder/QuotientRemainder (poly_divrem_modulus). VERIFIED.
- [x] A7: Discriminant degree-0/1 → 1 (poly.c). VERIFIED.
- [x] Tests: test_polymod.c (A6) + test_poly.c test_polynomial_modulus (A2/A3), test_discriminant (A7). PASS.

## Phase 3 — FLINT nmod_poly fast path (section C) — DONE
- [x] nmod_poly bridge (flint_bridge.c, #ifdef USE_FLINT): flint_nmod_poly_divrem, flint_nmod_poly_xgcd; word-size prime; classical fallback (poly_divrem_mod_classical).
- [x] Wired into Quotient/Remainder/QR (poly_divrem_modulus) + ExtendedGCD. GCD left as-is (out of A-list, variadic no-var API).
- [ ] check-c99 portability pass (deferred to end-of-phase gate)

## Phase 4 — remaining A-series (A4, A5, A8, A9, A10) — DONE (all verified + tests pass)
- [x] A4: PadRight 1-arg ragged — pr_is_atomic guards pr_scan_dims + pr_build (src/list/pad.c). Non-List-head padding preserved.
- [x] A5: Lookup Key unwrap — builtin_lookup (src/assoc.c). Threads over assoc-lists; Missing on absent.
- [x] A8: Transpose top-two-levels 1-arg (get_top2_dims/transpose_top2, src/list/transpose.c). n-arg unchanged.
- [x] A9: OptionValue under OptionsPattern[other] ($OptionsPatternHead$ in match.c + symtab.c). outer[]->111.
- [x] A10: Function slot-form closure — function_is_slot_form gates substitute_slots (src/purefunc.c).
- [x] Tests: padright/association/list/options/purefunc suites all PASS.

## Phase 5 — behavioural (B3, B2, B1)
- [x] B3: Print pipe-flush — fflush after each Print (src/print.c builtin_print). VERIFIED.
- [x] B2: 1/(4 Sqrt[2]) denominator grouping + numeric-first ordering (src/print.c Times). VERIFIED std+InputForm; no regressions.
- [ ] B1: ToNumberField canonical primitive element (assessing feasibility/scope)

## Phase 6 — remove .m workarounds + verify — DONE
- [x] ReduceModP/PolyQuoModP/PolyRemModP/PolyExtGCDModP delegate to fixed builtins (A2/A3/A6 + FLINT hot path)
- [x] PadRows -> PadRight[rows] (A4); index loops -> Do[..,{p,ps}] (A1); Lookup Key (A5); Options[iPIM] removed (A9)
- [x] Left as-is: Thread (A8, works correctly) and B1 theta-tolerant adaptation (needs B1 fix)
- [x] test_parallelmixedtower.c PASS incl. test_nontorsion_divisor_certificate; Coth[x]/(1+Sech[x]^5)^(3/2) now certifies non-elementary

## Verification gate
- [x] Targeted unit suites: iter, poly, polymod, list, association, options, purefunc, expr, core_algebra, packed_list — all PASS
- [x] check-c99 clean; check-packed-aware clean
- [ ] Full unit suite (running; dsolve_corpus last — unrelated to changes)
- [ ] valgrind on new-code paths (polymod/iter)
- [ ] B1: DEFER (deep FLINT qqbar canonicalization; flag to user)
- [ ] version.h bump + changelog (changelog DONE) + commit/tag (user decision)

## Changelog
docs/spec/changelog/2026-09-21.md — divergence-fixes section added.

## Review
A1–A10 correctness + FLINT nmod_poly fast path + B3 landed and pushed (v0.170,
commit 3c12c301). Broad backstop suite (476 pass / 18 fail) then classified the
18 against a baseline build of the parent (f9670b19):
- 8 fails were B2 printing-rot (my correct new output; tests encoded the old
  1/4/Sqrt[2] form). B2 reverted (v0.171, 57ea1f8e) and deferred as a dedicated
  test-expectation sweep, alongside B1; those 8 suites now pass.
- 5 fails are PRE-EXISTING (fail identically on baseline f9670b19), NOT caused by
  this work: risch_rde_tower (Integrate[E^(Log[x]^2)] FreeQ=True), intrischnorman
  (1/Log[x] -> ExpIntegralEi vs LogIntegral), moebiusmu (flaky: ECM
  non-determinism on 10^50+1), primenu (baseline also 7), qrdecomposition_machine
  (segfault on baseline too — gcc-16 class).
- 4 slow corpus/dsolve suites (crc_corpus, intrat_corpus, dsolve, dsolve_stress):
  verifying on v0.171 (were 90s-timeout false-fails and/or B2-rot now reverted).

Lesson: a print-formatting change has a broad blast radius via test expected-output;
run the full suite (or at least all print-touching suites) BEFORE pushing such a
change to main. Pushed before the backstop completed — corrected with the revert.
