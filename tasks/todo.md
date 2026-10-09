# Task: Simplify review — faster, more robust, less prone to hanging

Plan: `~/.claude/plans/i-would-like-to-binary-pie.md`. Sequencing
1(termination) → 2(leak) → 3(score memo) → 4(0-backstop) → 5(gates).

## Progress

### Phase 1 — Guarantee termination  ✅ DONE (all tests green)
- [x] 1.1 `pseudo_rem_standard` (`poly.c:4538`) monotonicity guard (`prev_degR`
      + `stalled` flag skips the lc(B) padding on a stall break).
- [x] 1.2 Whole-call `TimeConstraint` deadline (`g_simp_call_deadline` in
      `simp_search.c`): armed once in `builtin_simplify`, honoured in
      `simp_dispatch`, `simp_bottomup`, the three `simp_pipeline_*`, the seed
      phase (mid-phase check before the radical cluster) and the round loop
      (`simp_past_deadline`). Scalar `t` now bounds the WHOLE call (closes the
      SHAPE_RATIONAL / bottom-up bypass); `{tLoc,tTot}` honoured
      (`simp_parse_total_budget`).
- [x] 1.3 Monomial-reduction loop (`poly.c:~5372`) size+iteration guard.
- [x] Regression: 14 Simplify suites + 8 poly/factor suites all green.

### Phase 3 — Score memoization  ✅ DONE (behaviour-identical; 22 suites green)
- [x] `update_best_scored(best,score,c,s)` + `consider_candidate(...)` helpers in
      `simp_search.c`; each round-loop candidate scored ONCE (was 2–4×). Folded
      5 near-identical extra-candidate blocks (TrigRoundtrip / PythagSquareComplete
      / PythagReduce / HalfAngle / Radicals) into one helper. Main `r`: scored once,
      reused for update_best + propagation gate + TrigExpand bound.
- [x] A/B speedup vs clean `main`: pending (baseline building).

### Phase 4 — Literal-0 backstop  ❌ DROPPED (unsound) — reverted
- Attempted: on built-in `best == 0`, re-verify with `zero_test_decide_assuming`
  and override on ZERO_TEST_FALSE. REVERTED: `zero_test_decide` has false
  negatives (FALSE = "proved OR strongly believed nonzero"); it returned FALSE on
  the genuine identity `(2+Sqrt[5])^(1/3)+(2-Sqrt[5])^(1/3)-1 == 0`, corrupting a
  correct 0 (regressed `radical_simplify_tests`). User confirmed the rule:
  **Simplify must not use PossibleZeroQ** ([[feedback_simplify_must_not_use_possiblezeroq]]).
  The FactorTerms→0 class is already fixed at source (ft_content_wrt_set unit
  content), which is the correct place.

### Phase 5 — Regression gate  ◑ PARTIAL
- [x] Extended `test_simplify_hang.c`: `{tLoc,tTot}` whole-call cap coverage,
      rational-shape dispatch-path option acceptance, stale banner fixed.
- [ ] `tests/bench_simplify.c` perf/leak gate — DEFERRED with the leak fix (its
      leak tripwire can't go green until the leak is fixed).

### Phase 2 — Leak fix  ⏸ DEFERRED (per user: finish 3–5 first)
- Reproduced under valgrind: `Simplify[2 Cos[x - y/2] Sin[3 y/2]]` leaks
  ~11 definitely-lost blocks/call (N=1→431, N=40→860). Localised to the
  round-loop transform site `simp_search.c:1922` (the leaked block is part of
  `r`, a TrigToExp output tree; grows 8→367 blocks N=1→40).
- RULED OUT: the FactorMemo path (forcing memo off in `trig_to_exp` left the
  per-call leak unchanged); the three chained passes take `const Expr*` and
  free their own results; the CandSet lifecycle (`cs_add_or_free`/`cs_free`)
  and the `seeds=next` transition are balanced; standalone
  `Expand[TrigToExp[...]]` does NOT leak (no per-call growth).
- STILL UNPINNED: `r`'s tree is lost only in the full Simplify context; root
  cause resists static analysis — needs multi-cycle runtime bisection
  (selectively disabling transforms/passes under valgrind).

### NEW FINDING (out of plan scope) — multi-trig evaluator hang
- `Simplify[Sum[Sin[k x]^4 Cos[k x]^2 + Tan[k x]^3, {k,1,n}]]` hangs even for
  small n. `sample` shows the cost is in the EVALUATOR, not Simplify:
  `evaluate_step → builtin_times → trig_canon_groups` (`src/trig_canon.c`).
  No Simplify-layer `TimeConstraint` can interrupt a kernel inside `evaluate()`.
  Pre-existing (my changes don't touch eval/times/trig_canon). This is a real
  user-visible "Simplify hangs", rooted in `Times` trig canonicalisation.

### Phases 3–5 — not started
- [ ] 3 score memoization (pure perf; low risk)
- [ ] 4 literal-0 backstop (robustness)
- [ ] 5 regression gate (perf + termination)

## Review (2026-10-09)

**Outcome vs the three goals:**
- *Less prone to hanging / termination:* real, verified wins — the last unguarded
  pseudo-division loop is closed, and a user-set `TimeConstraint` now bounds the
  whole call on every path (was ignored on rational-shape / bottom-up inputs).
- *More robust:* the literal-0 backstop was dropped as unsound; the net robustness
  gain is the termination guarantees + keeping Simplify free of PossibleZeroQ.
- *Faster:* **no measurable wall-clock speedup found at the search layer.** A/B
  against clean `main` on a scoring-heavy custom-`ComplexityFunction` workload
  (x3000): baseline ~10.3 s vs after ~10.3 s. Simplify's cost is dominated by the
  transform kernels and the evaluator, not scoring/classification. Phase 3 is kept
  as a clean behavior-identical refactor (removes redundant 3×→1× scoring, dedups
  5 blocks), not a perf win.

**Where the real speed/hang levers are (for a future pass):** the evaluator
`trig_canon_groups` (`builtin_times`, multi-distinct-argument trig — a genuine
hang), and the heavy kernels (multivariate `Together`/`Factor`). Both are outside
`src/simp/` with broader blast radius.

**v0.326** behaviour change = the whole-call `TimeConstraint`. Changed tracked
files: `src/poly/poly.c`, `src/simp/{simp_search.c,simp_builtins.c,simp_bottomup.c,
simp_internal.h}`, `src/version.h`, `tests/test_simplify_hang.c`,
`docs/spec/builtins/simplification.md`, `docs/spec/changelog/2026-10-05.md`.

## Notes
- Committed as `v0.326` and tagged `v0.326`.
- Deferred, tracked: the leak (Phase 2), the `trig_canon_groups` evaluator hang,
  and `tests/bench_simplify.c` (its leak tripwire needs the leak fixed first).
