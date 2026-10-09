# PossibleZeroQ backlog #4 + #5 (v0.323)

Finish the two OPEN items in `POSSIBLE_ZEROQ_IMPROVEMENTS.md`. They ship
together (fixing #4 re-exposes #5). Plan: `~/.claude/plans/what-s-left-to-do-radiant-journal.md`.

## Fix #5 — cancellation-aware screen FALSE (lands first)
- [ ] Add `#define ZT_HUGE_SCALE_BITS 40` near the `ZT_*` tunables
- [ ] Add `huge_scale_confirms_nonzero(e, mag, scale)` helper (one MPFR rung, persist>=half)
- [ ] Gate the fast-FALSE in `screen_point` (:1171) — zero-ish branch returns TRUE
- [ ] Gate the fast-FALSE in `decide_numeric` rung-0 (:1102) — zero-ish continues ladder honored

## Fix #4 — re-draw a numeric-function-head decline (lands second)
- [ ] Add `residue_is_numeric_decline(z)` predicate (ATTR_NUMERICFUNCTION walk)
- [ ] Broaden the `!is_pure_numeric` branch in `evaluate_rung` (:1048) to set re-draw
- [ ] (optional) rename `out_overflow` -> `out_redrawable` + update comments
- [ ] Confirm Gamma/Erf/Erfi/Fresnel carry ATTR_NUMERICFUNCTION

## Tests
- [ ] New Group 19: Gamma[x^2]+1 -> False; Erf sibling; determinism x4; UndefinedZQHead trip-wire
- [ ] Robustify `test_battery_weierstrass_cosh_product_roundtrip` across seeds
- [ ] Add genuine huge non-zeros: Tanh[20]+2^60 -> False; 10^20+1 -> False
- [ ] Fix Group 10 header overclaim

## A/B regression (restore baseline zero_test.c, NOT git stash)
- [ ] zero_test / possiblezeroq_* / trigexp_zero suites identical except intended flips
- [ ] dsolve-verify slice, comparisons, simplify/solve/limit unchanged

## Housekeeping
- [ ] make -j + make check-c99 clean
- [ ] REPL spot-checks (Gamma[x^2]+1, Weierstrass, Tanh[20]+2^60, Pythagorean identity)
- [ ] POSSIBLE_ZEROQ_IMPROVEMENTS.md: mark #4/#5 RESOLVED
- [ ] src/version.h -> v0.323
- [ ] docs/spec/changelog/2026-10-05.md note
- [ ] docs/spec/builtins zero-test page note
- [ ] commit (; v0.323) + tag v0.323 + push --follow-tags

## Review (v0.323 — shipped #4 only; #5 found non-reproducing)

**Outcome.** The backlog #4/#5 example repros no longer trigger on HEAD (Gamma was
independently fixed to return finite extended-range values; #5's deep-cancellation
FALSE is masked by poles and machine accuracy). BUT the #4 *mechanism* is alive via
`PolyLog`: `N[PolyLog[2, z]]` declines for `|z| >~ 300`, so
`PossibleZeroQ[PolyLog[2, x^2] + 1]` returned a silent wrong `True`. Fixed it.

**Shipped (src/zero_test.c, +58/-9).**
- `residue_is_numeric_decline(z)` — a residue built only from `ATTR_NUMERICFUNCTION`
  heads on numeric args is a magnitude decline, not a genuine symbolic residue.
- `evaluate_rung` flags such a residue re-drawable (reuses the `out_overflow`
  channel); `sz_trial_shelled`'s shell ladder shrinks the sample to a resolvable
  point → decides FALSE. Undefined heads (no NumericFunction) still abort → True.
- Comments at the predicate / `evaluate_rung` / `sz_trial_shelled` / Phase-A screen
  updated to the PolyLog live case.

**Dropped.** The #5 huge-scale cancellation guard: prototyped, measured against a
pre-change binary, changed NO verdict (unverifiable on current inputs); the user
chose "fix the real bug" over "ship the hardening." Documented in the backlog.

**Verified.** zero_test_tests (incl. new Group 19), possiblezeroq_* (3), trigexp_zero,
comparisons, simplify, solve, limit, refine, interp — all PASS. dsolve_tests fails
PRE-EXISTING (`SpecialFunctionForm` assertion, identical on baseline A/B; unrelated).
make check-c99 + check-messages clean. A/B: PolyLog cases flip True→False, identity
and trip-wire unchanged.

**Housekeeping done.** version.h → 0.323; POSSIBLE_ZEROQ_IMPROVEMENTS.md #4 RESOLVED
/ #5 non-reproducing note; docs/spec/changelog/2026-10-05.md; expression-information.md
re-draw paragraph. NOT yet committed/tagged (awaiting go-ahead).

**Process lesson.** An A/B that rebuilds two source variants into the SAME binary
path (`./Mathilda`) silently compares a binary against itself if you forget which
build is currently there — it produced a wrong "dead code" conclusion mid-session
until caught. Build the two variants to DISTINCT paths.
