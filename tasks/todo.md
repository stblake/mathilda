# NSum engine deficiency fixes — DONE 2026-06-14

Plan: `/Users/user/.claude/plans/cuddly-napping-orbit.md`. All changes in
`src/numerical_calculus/nsum.c` + tests + docs.

## D1 — peaked / late-settling summands ✅
- [x] Far-tail magnitude ladder (run only when head isn't already monotone)
- [x] `ns_diverges` → sustained far-tail growth (`rising_far`)
- [x] decay-onset `settle` → EM adaptive base `N=imin+max(NSumTerms,settle)·di`
- [x] Verified: 1/(1+(k-20)^2)→3.10462; 2^k still CI; harmonic not flagged

## D2 — accuracy ceiling ✅
- [x] dequad reltol/levels scaled to target precision (just above the floor)
- [x] Hybrid corrections: symbolic D primary, numerical contour DFT on balloon
- [x] Guard digits (2·target+32 bits) vs near-1 float-sample cancellation
- [x] Oscillatory (sign-alternating) extensions stay symbolic → NULL → Wynn
- [x] Verified: Log[1+1/n^2] & EM-gamma @WP35 → ~33–34 digits (were 8–12)

## D3 — nested summands ✅ (routing) / documented residual
- [x] Removed wrong-premise black-box EM-exclusion; EM valid for independent inner
- [x] Verified: multidim {k,1,∞}→1.14434, {k,1,n}→0.770188
- [ ] (residual, documented) infinite-outer multidim PRODUCT mixed series ≈0.607 vs 0.564

## Wrap-up ✅
- [x] `test_peaked_summand` + `test_composite_precision` added to test_nsum.c
- [x] nprod WP composite tolerance tightened 1e-11 → 1e-25
- [x] nsum/nprod/nseries/nlimit tests pass; suite under the 60s alarm
- [x] valgrind matches pre-change baseline (no new leaks)
- [x] docs/spec calculus.md + changelog 2026-06-08.md + memory updated

## Review
The plan's "robust full fix" (replace symbolic EM derivatives with contour) was
too aggressive on its own: contour is fragile for geometric/oscillatory summands
and the black-box exclusion broke valid nested EM. The shipped design is a
**hybrid** — symbolic D primary (byte-identical to the original for simple
summands, so no regressions), contour only when symbolic balloons (composite
summands), never for oscillatory extensions — plus WP-scaled integral tolerance
and guard digits. This fixes D1 and D2 fully and corrects D3's routing, with one
genuinely-hard nested mixed-series case left documented. Key trap: NSum's shared
evaluator has a pre-existing per-eval GMP-rational leak, so any added per-call
summand evaluation amplifies it — the far-tail ladder is skipped on monotone
heads and oscillation is read from existing head terms to keep evals flat.
