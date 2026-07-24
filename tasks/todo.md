# Gruntz `Limit` Phase 3 — deep log-tower cancellation (thesis 8.19)

Per GRUNTZ_STATE.md "Known gaps"; the chosen next Gruntz gap.

- [x] Diagnose: `Method->"Gruntz"` abstains on 8.19; failing atom is the
      leading-order cancellation `Log[x+Log[x]]-Log[x]`.
- [x] Fix 1 (`expand_logs`): factor the `w`-pole out of each `Log`
      (`Log[G]=k Log[w]+Log[H]`, `H=Expand[G w^{-k}]`) so `Series` never splits
      it into the branch-contaminated `Log[-a]+Log[-1/a]` (`=2 Pi I`).
- [x] Fix 2 (`series_leadterm`): run `Series` in the positive log-scale
      `P=-Log[w]` (substitute `lw_sym->-P`, restore `P->-Log[w]`) so
      `Log[-Log[w]]` stays real; drop the now-unnecessary `Log[lw_sym]` bail.
- [x] Fix 3 (`series_leadterm`): accept a `w`-free level as the sub-scale
      coefficient (`limitinf` recurses) instead of abstaining.
- [x] Tests: new `test_log_tower` (atoms + full 8.19 `->1`); 8.19 removed from
      `test_honest_abstentions`; `gruntz_tests` green.
- [x] Regression: `limit`/`limit_assumptions`/`nlimit`/`nseries` 0 FAIL;
      `series`/`residue` fail-modes IDENTICAL to pre-change (ArcCoth family +
      NResidue convergence warnings, both pre-existing).
- [x] valgrind: definitely-lost 13,440B/420blk == documented macOS baseline; no
      new gruntz frames in any leak stack.
- [x] Docs: `GRUNTZ_STATE.md` (Phase 3 + gap update), `calculus.md` (Gruntz
      coverage), weekly changelog `2026-07-20.md`.

# Gruntz `Limit` Phase 3b — fix the Automatic/Series cascade hang on 8.19

Follow-up to Phase 3 (user: "triage and fix the asymptotic series hangs").

- [x] Triage: plain `Limit[8.19]` (and `Method->"Series"/"Asymptotic"`) hung.
      Instrumented `limit.c`'s `TRY` cascade -> the recursion reaches the Gruntz
      layer at depth ~23 on a `x->1/x`-substituted deep sub-limit and hangs
      there (not in the asymptotic layers themselves).
- [x] Root cause: `series_leadterm` fed `Series` a complex-contaminated `f4`
      (`P + I Pi ...`). A 3+-level tower re-introduces `Log[-Log[w]] ->
      Log[Log[w]] + I Pi` below the single positive-scale substitution; Series
      loops on it.
- [x] Fix: `contains_complex_head` guard — reject a frozen-scale Series input
      carrying `I` (mrv engine is real-valued, so it is always a branch
      artefact). `series_leadterm` abstains; `expand_logs` leaves such a log
      unexpanded. Fast abstention on the sub-limit; top-level Gruntz then closes
      8.19 -> 1.
- [x] Fixed a double-free introduced mid-work (`mk_log(G2)` consumes G2; the
      `ser` complex-guard must not also free G2).
- [x] Tests: `test_log_tower_no_hang` (Automatic 8.19 -> 1; Series/Asymptotic
      terminate). Un-hangs `test_series_infinity_no_inv_var_leak`, so
      `series_tests` now completes.
- [x] Stress: 21-case log-tower battery (`scratchpad/stress.py`) — 0 hangs, 0
      wrong values (one deep 3-tower is a bounded ~15s honest abstention).
- [x] Regression: gruntz/limit/limit_assumptions/nlimit/nseries green; series
      failures are a SUPERSET of committed (proven: the 3 extra D/Integrate ones
      are pre-existing, previously hidden behind the hang). valgrind == baseline.

## Review (Phase 3b)

**Net effect.** Plain `Limit[8.19]` now resolves to 1 instead of hanging, and
the fix un-hangs `series_tests` (which contained a plain-`Limit[8.19]` test).
The key insight: the mrv engine is real-valued, so any `Complex[]` surviving in
a frozen-scale `Series` input is a branch artefact — rejecting it both prevents
the hang and preserves "never wrong" (the engine abstains rather than returning
a complex/garbage value). Verified there is no regression: the 3 newly-visible
`series_tests` failures are pre-existing `SeriesData` D/Integrate normalisation
gaps that the committed suite never reached because it hung first.

---

# Gruntz `Limit` Phase 3 — deep log-tower cancellation (thesis 8.19), review

**What shipped.** Three surgical fixes confined to `src/calculus/gruntz.c`
(`expand_logs`, `series_leadterm`) make the `Method->"Gruntz"` engine resolve
thesis 8.19 and the class of pure log-tower limits with leading-order
cancellation. Root cause was two ways `Series` mishandles the frozen
`lw_sym = Log[w]` (whose sign it doesn't know): it splits a `w`-pole log into a
spurious `2 Pi I` constant, and it mis-branches `Log[-lw_sym]` (a real,
positive-scale log) to `Log[lw_sym] + I Pi`. Factoring the pole out first and
expanding in the positive scale `-Log[w]` both address that at the source; the
third fix lets a fully-collapsed level recurse as its sub-scale limit.

**Scope honesty.** Only the explicit `Method->"Gruntz"` reaches the fix. Plain
`Limit` (Automatic) still hangs on 8.19 inside `limit.c`'s Asymptotic/Series
cascade layers, which run before the Gruntz layer and never call `gruntz_limit`
— a pre-existing issue (byte-identical before/after this change), left as a
separate follow-up.

---

# FLINT-backed multivariate zero test (`flint_mpoly_is_zero`)

Implements FLINT_ZERO_TEST_PLAN.md.

## Phase 1 — Primitive
- [x] `flint_bridge.h`: declare `int flint_mpoly_is_zero(const Expr* e);`
- [x] `flint_bridge.c` (USE_FLINT): kernel-aware `to_mpoly_kernels` + `flint_mpoly_is_zero`
      - operate on `expr_expand(e)` (mirrors classical; captures Tan·Cos→Sin reductions)
      - generator set = poly.c `collect_variables` (symbols + opaque kernels)
      - match-generator-first, then Plus/Times/Power(nonneg int)/Rational/Integer/Bigint
      - decline (-1) on anything FLINT can't model (Complex, Real, Pi, frac power, …)
- [x] `flint_bridge.c` (!USE_FLINT): stub returning -1

## Phase 2 — Wiring + differential harness
- [x] `poly.c`: `is_zero_poly` dispatches to FLINT path, classical `is_zero_poly_depth` fallback
- [x] env-gated (`MATHILDA_ZEROTEST_DIFF`) differential harness: run both, abort on disagreement

## Phase 3 — Tests
- [x] extend `tests/test_zero_test.c`: kernel-generator zeros/non-zeros, decline cases
- [x] random-polynomial fuzzer under differential mode (400 iters, 0 disagreements)

## Phase 4 — Verify
- [x] Build clean (make + cmake tests)
- [x] Differential run: 1 disagreement found (Tan·Cos) → fixed (expand-first); zero_test+simplify clean
- [x] Regression suites (20 core + 9 integrate/corpus): 0 disagreements everywhere;
      simplify/ratcanon_spec/intrat/series/intrat_corpus/goursat/cherry_ei fail-modes
      IDENTICAL to pristine (all pre-existing, verified by stash+rebuild)
- [x] valgrind clean — definitely-lost 13,440B/420blk == documented macOS baseline; no
      leak stack in new code (flint_mpoly_is_zero/to_mpoly_kernels/collect_variables)

## Phase 5 — Gruntz Gamma acceptance (stretch) — DONE
- [x] re-enabled Gamma isolation in gruntz.c (Gamma[g]=Exp[LogGamma[g]] branch)
- [x] Stirling ratio -> 1 (0.06s), thesis 5.5 -> Infinity (0.81s) NOW RESOLVE
      (promoted to test_thesis_gamma); 8.31 still abstains (residual Series cost)
- [x] BENCHMARK: 8.31 abstention 173.9s (classical zero) -> 14.9s (FLINT) = ~11.7x
      — directly confirms is_zero_poly was the mrv bottleneck that forced Gamma out
- [x] gruntz_tests green

## Phase 6 — Docs — DONE
- [x] docs/spec/builtins/expression-information.md — PossibleZeroQ Stage-2 FLINT note
- [x] flint_bridge.h contract comment + weekly changelog (2026-07-20.md)
- [x] memory note updated (project_flint_zero_test DONE) + MEMORY.md pointer
- [x] GRUNTZ_STATE.md — Gamma re-enabled section

## Review

**What shipped.** `flint_mpoly_is_zero` (`src/poly/flint_bridge.{c,h}`) — a
FLINT `fmpq_mpoly`-backed accelerator for `is_zero_poly`, wired into
`src/poly/poly.c` with the classical `is_zero_poly_depth` as the decline
fallback and an env-gated (`MATHILDA_ZEROTEST_DIFF`) differential harness. Plus
the Gruntz `Gamma` isolation it unblocks (`src/calculus/gruntz.c`).

**Key design decision.** The accelerator operates on `expr_expand(e)`, not the
raw tree — the classical predicate does too, and expansion performs kernel
reductions (`Tan[k]·Cos[k]→Sin[k]`) that the raw generators hide. The
differential harness *caught* this as a real disagreement during bring-up; the
expand-first design makes the verdict bit-for-bit classical while still removing
the multiplicative re-expansion fan-out.

**Verification.** 0 differential disagreements across 29 suites + a 400-iter
fuzzer; all non-green suites verified fail-identical to pristine (stash+rebuild);
valgrind at macOS baseline; clean strict-C99 build. Gamma benchmark
173.9s→14.9s confirms the diagnosed bottleneck.

**Scope honesty.** The `p - Expand[p]` micro-bench shows no speedup (its cost is
the outer Expand). The real win is the deep-cancellation multi-generator mrv
towers — exactly the Gamma cases. 8.31 still abstains (residual `Series` depth,
out of scope for this zero-test project).
