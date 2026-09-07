# M21 — §2.2.1 corpus (Problems 1–100) with initial-condition coverage

Goal: add Nasser Abbasi §2.2.1 (100 elementary ODEs, 63 IVPs) to the DE corpus with real
initial-condition verification, then chase full 100/100 coverage by fixing DSolve's IVP
constant-fitting for transcendental-inverse / Riccati forms. 0 FAIL invariant.

## Stage 1 — Converter: emit + normalize ICs (`tools/latex_ode_to_mathilda.py`)
- [ ] Emit ICs as a DSolve-native equation list `{ode, ic1, …}` (stop discarding `conds`)
- [ ] Symbolic ICs `y(a)=b` → `y[a]==b`; derivative ICs `x'(0)=10` → `x'[0]==10`
- [ ] indVar fallback for swapped-variable rows #98–100 (x=x(y))
- [ ] Guard `is_condition_row` against the main primed ODE
- [ ] Test-convert all 100 → 100 non-placeholder rows, 0 "NO ODE ROW"

## Stage 2 — Prelude IC verification (`DSolve_test_status/dsolve_corpus_prelude.m`)
- [ ] `dsolveCheckCode` passes list `eqn` to DSolve; system guard still keys on `ListQ[fn]`
- [ ] `dsBranchVerdict` verifies every member (ODE residual + each IC), require all OK
- [ ] Symbolic-IC record verifies (a,b instantiated by dsFreeParams)

## Stage 3 — Generate corpus, register ctest, measure baseline
- [ ] Generate `DSolve_test_status/DE_examples_221.m` (--label 2.2.1); spot-check 10
- [ ] Register `dsolve_corpus_2_2_1_tests` in `tests/CMakeLists.txt`
- [ ] Full §2.2.1 TSV run → `reports/2.2.1.md` (pre-fix baseline + ranked failures)

## Stage 4 — Solver fixes (chase full coverage)   [pre-fix §2.2.1: 86/100]
- [x] Swapped-variable `A/x'==B` #98/99/100 — NthAlgebraic clears top-deriv denominator
      (Numerator[Together] + recurse). Worked around Exponent bug (see lessons) via FreeQ.
- [x] Transcendental-inverse IC fit — `dsolve_fit_constants` uses SCALAR Solve for the
      single-condition/single-constant fit (list-form Solve skips inverse-function
      solving). Fixes #33/#34/#61/#29/#30, and #40 (Airy Möbius fit).
- [x] `ConditionalExpression` principal-branch collapse (multivalued Tan inverse) — fixes
      #60 (y'=1+y²,y0=0 → Tan[x]).
- [ ] #67 malformed separable body (spurious 2πI C[1] + tangled Log) — investigate
- [ ] #35 `y'=Log[1+y²]` hangs → bounded UNEVAL (non-elementary; matches Maple/Mma) — OK
- [ ] #47/#48 messy Root/slow — check if they PASS or stay bounded UNEVAL
- [ ] No regression: `dsolve_corpus_2_1_2_tests` @572, all dsolve_* suites, make check-c99

## Follow-up bugs found (separate from M21)
- `Exponent[expr, v]` returns 0 whenever expr contains ANY funcapp (Sin[y], x[y]).
  Worked around in NthAlgebraic with FreeQ. Real bug in src/poly/exponent.c — fix later.
- Corpus prelude `dsFreeParams` had `Heads->True` → collected Plus/Times/Tan as "params"
  and substituted numbers for them → every residual non-numericizable → vacuous UNK.
  FIXED (removed Heads->True); this makes §2.1.2 numeric verify real for the first time.

## Stage 5 — Re-measure, gate, document
- [ ] Set `dsolve_corpus_2_2_1_tests` baseline to achieved high-water
- [ ] STATUS.md §2.2.1 block + wave-history line; reports/2.2.1.md
- [ ] DSOLVE_PLAN.md M21 entry; changelog 2026-09-07.md; README.md recipe update

## Review

**Result: §2.2.1 86 → 96 / 100 scalar PASS, 0 FAIL, 0 regression.** All with real
initial-condition verification (ODE residual + every IC).

Done:
- Stage 1–3: converter emits ICs (DSolve-native `{ode, ic...}` list), swapped-var + symbolic-IC
  handling, section-agnostic table/column detection. Prelude verifies every member; a
  leaked-constant IVP → UNEVAL. `DE_examples_221.m` (100 recs, 63 IVPs) + ctest #241 (baseline 4)
  + `reports/2.2.1.*` + STATUS.md block.
- Stage 4 fixes: `NthAlgebraic` denominator-clearing (98/99/100); scalar-Solve IC fit
  (29/30/33/34/40/61); `ConditionalExpression` principal branch (60). #36 regression caught & fixed
  (scalar-empty ≠ no-solution at singular fit point).
- Verifier bug fixed: prelude `dsFreeParams` `Heads->True` → vacuous UNK. Now real.

Verification:
- dsolve_tests + dsolve_stress_tests: **all pass** (no regression from C changes).
- make check-c99: **clean**.
- §2.1.2 regression (final binary): **436/1000 PASS (+4 vs 432), 564 non-PASS ≤ 572 gate,
  0 FAIL.** Shared fixes helped 4 §2.1.2 cases; verifier fix surfaced no new FAIL.
- `ctest -R dsolve_corpus_2_2_1_tests`: **Passed** (baseline 4).

Residue (4, bounded UNEVAL, 0 wrong answers): 35 (non-elem + missed equilibrium), 47 (deg-12
Root), 48 (slow ArcSin), 67 (Solve C[1] branch/integration-constant collision).

Follow-up bugs filed (memory + STATUS): `Exponent[expr,v]`=0 with any funcapp; Solve
generated-constant collision.
