# Gruntz stress-test generalization + advanced tutorial

## Goal
- For each Gruntz thesis example (8.1–8.37 + intro/Ch.2 worked examples), produce a
  family of 10+ increasingly-difficult *generalized* limits with hand-derived exact
  answers.
- Verify each against Mathilda `Limit[..., Method->"Gruntz"]`; keep only those whose
  output is mathematically correct. Assemble passers into `tests/test_gruntz_stress.c`.
- Turn the verified families into an advanced tutorial page under the MkDocs site
  tutorials section (calculus/limits).

## Plan
1. [x] Build Mathilda; NDJSON batch harness (scratchpad/glimit.py, per-case isolated).
2. [x] Generate candidate families (scratchpad/gen.py, gen2.py) with hand-derived answers.
3. [x] Run all candidates; keep verified passers. 141 verified across 11 families.
4. [x] Assemble tests/test_gruntz_stress.c; wire into tests/CMakeLists.txt; green (0 FAIL).
5. [x] Write advanced tutorial page; wire into .pages nav; cross-link from 08-calculus.
6. [x] Changelog + GRUNTZ_STATE update.

## Verification rule
Pin ONLY outputs Mathilda actually returns AND that I have hand-verified correct.
Abstentions/wrong answers are dropped (optionally pinned as honest-abstention docs).

## Review
- **Deliverables:** `tests/test_gruntz_stress.c` (new `gruntz_stress_tests` suite,
  141 verified generalizations + 4 honest abstentions, all green under Release);
  `site/docs/tutorials/computing-limits-gruntz.md` (advanced tutorial, 53
  transcripts verified exact vs. binary); nav + cross-link + changelog + state.
- **Families (base → axis):** A max-base 8.12; B cancellation 8.1 (→ -c); C
  exp-tower 8.5–8.8 (→ E^a); D trig-at-vanishing 8.21/8.22; E nested-log 8.19/8.20;
  F Hardy scale 8.9; G conjugate radicals 2.5; H finite-point power series 2.6/2.7;
  I special-fn singularities 8.23–8.34; K Max/Min 8.37; L digamma/log-Gamma growth.
- **Verification wins:** isolated harness exposed a false "hang" batching artifact;
  8.10-as-transcribed → Infinity (dropped, not the thesis's 1/3); Ei[2x]… = 1/2 not 1.
- **Honest gaps confirmed:** 3-level tower, Gamma/psi *difference* asymptotics,
  bare Gamma ratio — pinned as abstentions, never a wrong value.
- **No engine code changed** — tests + docs only; existing gruntz_tests still green.
