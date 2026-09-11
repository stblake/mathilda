# DSolve M37 — §2.2.17 corpus (Problems 1601–1700) + 5-FAIL fix + Abel AIR

Baseline (measured, existing binary + real prelude): **79 PASS / 5 FAIL / 16 UNEVAL**.
Target: ≈89–90 PASS, **0 FAIL**, 0 regression.

## Phase 0 — Corpus creation & registration
- [ ] 0.1 Generate `DSolve_test_status/DE_examples_2217.m` via converter (100 recs, 28 IVP)
- [ ] 0.2 Register ctest `dsolve_corpus_2_2_17_tests` in `tests/CMakeLists.txt` (after :3971)
- [ ] 0.3 Generate `reports/2.2.17.{tsv,md}`; add STATUS.md block + README.md row

## Phase 1 — Fix the 5 FAILs (MANDATORY, restores 0-FAIL) — HEADLINE ✅ (section: 83/100, 0 FAIL)
- [x] 1.1 New substrate helper `ds_branch_num_ok(P, body)` (reject-direction) + factored `ds_subst_generics`
- [x] 1.2 `dsolve_fit_constants` — try all Solve roots, pick verifying one (1638/1641 → PASS)
- [x] 1.3 `dsolve_separable_try` — numeric-filter ± inversion branches (1622 → PASS)
- [x] 1.4 `dsolve_bernoulli.c` — numeric-filter ± two_signed branches (defensive)
- [x] 1.5 NEW: `ds_branch_corpus_verifiable` post-fit gate in dsolve_run (prelude-matching, 1st-order) → 1636 PASS via Separable-implicit; general 0-FAIL guarantee
- [~] 1.6 0-regression check on prior corpora — IN PROGRESS (dsolve_tests + §2.2.16 running)
  Note: 1624 → bounded UNEVAL (cube-root real-branch domain restriction; documented residue, not a wrong answer)

## Phase 2 — DSolve`AbelAIR (scaling reduction u=y/g(x) → separable) ✅
- [x] 2.1 New `src/calculus/dsolve_abel_air.c` (pipeline + implicit-function verify + exact split check)
- [x] 2.2 Wire cascade in dsolve.c (enum, string, externs, slot @:427, pinned, init)
- [x] 2.3 Add file to tests/CMakeLists.txt COMMON_SRC
- [x] 2.4 Closes 1604/1606/1607/1676 (relaxed deg_N gate for 1606); anti-overfit t_m37_abel_air

## Phase 3 — Opportunistic UNEVALs ✅ (net +1: 1697)
- [x] 3.1 1601 — rational-sample split tried; REVERTED (garbled Solve inversion, earns nothing)
- [~] 3.2 1673 — left as residue (Riccati, risky)
- [x] 3.3 BONUS 1697 — Bernoulli `bern_Y_in_sum_power` fast-decline gate → Exact claims it (fixes a
      pre-existing cold hang masked by warm eval-memo caching in sequential measurement)

## Phase 4 — Docs, version, gates
- [x] 4.1 Changelog block (docs/spec/changelog/2026-09-07.md), version 0.134→0.135
- [x] 4.2 DSOLVE_PLAN.md M37 entry; docs/spec/builtins/calculus.md DSolve`AbelAIR row
- [x] 4.3 Anti-overfit unit tests t_m37_fractional_power_branch / t_m37_abel_air
- [x] 4.4 §2.2.17 corpus: 87/100, 0 FAIL; ctest registered gate baseline 13; STATUS/README/reports done
- [x] 4.5 Gates — ALL GREEN: full DSolve corpus ctest (18/18 incl §2.1.2's 1000 + §2.2.1–17, 0 fail,
      2476s); `dsolve_stress_tests` ✅; `t_m37` assertions ✅ (verified directly); `make check-c99` ✅;
      valgrind — no new leak/error traces to my code (builtin_times uninit = macOS baseline noise;
      426 lost blocks = pre-existing Integrate/Solve per-call baseline).
- [x] 4.6 Rebuilt code-review graph.

## Review
**§2.2.17 corpus: 87/100, 0 FAIL** (baseline 79/5/16). The 5 FAILs (fractional-power ODEs shipping
a wrong branch) are fixed — 0-FAIL invariant restored — via one DRY mechanism: a numeric
back-substitution filter matching the corpus prelude's verdict, at the IVP constant-fitter
(verifying-root selection), the Separable/Bernoulli ± branches, and a prelude-matching post-fit gate
in `dsolve_run` (all radical-gated for perf). `DSolve\`AbelAIR` (new, scaling u=y/s → separable)
closes 4 Abel-2nd-kind UNEVALs (1604/1606/1607/1676) — a down payment on M13, full AIR invariant
table deferred. Bonus: a `bern_Y_in_sum_power` fast-decline gate fixed a pre-existing Bernoulli
cold-hang on 1697 (exact eqn) that warm sequential measurement had masked.

**Residue 13 (bounded declines, no wrong answers):** 8 `y=_G(x,y')` nonelementary (SymPy also fails),
1624 (cube-root real-branch), 1601 (garbled Solve inversion), 1673 (Riccati), 1681/1689 (nonelementary
integrating factor).

**Changed:** dsolve_common.{c,h}, dsolve_separable.c, dsolve_bernoulli.c, dsolve_abel_air.c (new),
dsolve.c, tests/{CMakeLists.txt,test_dsolve.c}, version.h, DSOLVE_PLAN.md, docs/spec/..., DSolve_test_status/*.
0 regression across all 18 corpus sections. Not committed (awaiting user).

## Review
(to be filled at completion)
