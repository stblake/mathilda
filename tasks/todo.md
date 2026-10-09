# Task: PossibleZeroQ `Assumptions` stress suite → fix coupling inequalities (#6)

User supplied a 20-case `Assumptions` stress suite (all identically zero on
their region, expected `True`) and asked to test it and improve `PossibleZeroQ`.

## Plan
- [x] Run the 20-case suite against HEAD; probe numericalization of each head.
- [x] Map the assumption-aware path in `src/zero_test.c` (sampler, specs,
      coupling classification).
- [x] Diagnose the single failure and design a sound fix.
- [x] Implement region-conforming rejection sampling (coupling inequalities).
- [x] Add `test_pzq_coupling_relations`; re-run suite → 20/20.
- [x] Regression-run zero_test / stress / comparisons / refine / reduce /
      fullsimplify / simplify_hang.
- [x] Audits: `make check-c99`, `make check-messages`.
- [x] Docs: spec `expression-information.md`, changelog `2026-10-05.md`,
      `POSSIBLE_ZEROQ_IMPROVEMENTS.md` #6; version bump v0.325; memory update;
      graph rebuild.
- [ ] Commit + tag `v0.325` + push — deferred to the user (not yet requested).

## Diagnosis (confirmed empirically against HEAD v0.324)
19/20 passed. The one failure: `PossibleZeroQ[Max[x, y] - x, Assumptions ->
x > y]` → `False` (should be `True`); `x >= y` variant identical; deterministic;
`N[Max[5,3]-5] = 0.0` so NOT a numericalization gap. Root cause: the
Schwartz–Zippel sampler draws each symbol independently from a per-symbol
`SampleSpec`, which cannot encode a fact coupling two symbols (`x > y`), and
`fact_keeps_false_sound` kept `False` for *every* inequality — unsound for a
**piecewise** head (`Max` is a different function off the region), so off-region
`x < y` draws reported a non-identity.

## Fix (src/zero_test.c) — region-conforming rejection sampling
A coupling dense relation (`>`,`>=`,`<`,`<=`,`Inequality`,`Unequal`) over ≥2 of
the expression's free symbols, groundable by the draw, is collected into
`conform_facts`. `sz_trial` re-draws any assignment outside the feasible region
(`assignment_conforms` → substitute + `bound_to_double` + compare on doubles; no
`Simplify`, no evaluator re-entry) up to `ZT_CONFORM_MAX_TRIES` (128); exhaustion
skips the point, and an all-skip screen returns honest `Unknown`. Equality
couplings (measure zero) keep the existing `False → Unknown` downgrade. New
helpers: `is_dense_relation_head`, `relation_true_binary`, `concrete_fact_holds`,
`assignment_conforms`, `collect_relation_vars`. The misleading "full-measure"
comment in `fact_keeps_false_sound` was corrected.

## Review (results)
- **Fix:** suite now 20/20; #15 `True` ×3 fresh procs (no flakiness).
- **Soundness (no over-downgrade):** `Max[x,y]-y` under `x>y` → `False`;
  `x-y` under `x>y` → `False`; `Refine[x==y, x>y]` → `False`; coupling equality
  `a-b` under `a==b` → `True` (downgrade preserved).
- **Byte-for-byte guarantee:** new path engages only when a ≥2-symbol coupling
  inequality is present ⇒ every non-coupling draw stream unchanged.
- **Regressions:** zero_test_tests, possiblezeroq_assumptions_tests (+ new
  group), possiblezeroq_stress_tests (104), comparisons_tests, refine_tests,
  reduce_tests, fullsimplify_tests, simplify_hang_tests — all pass. check-c99
  and check-messages clean.
