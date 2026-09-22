# TODO: DSolve M56 — 2nd-order integrating factors, Stage 2 μ(x,y′)

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-soft-turing.md`
Paper: Cheb-Terrab & Roche 1999, §2.2 (Lemma 3 Cases A/C/D + Lemma 2 μ̃).

Unblock M18 Stage 2 via path (b): emit the verified first integral
`R(x,y[x],y'[x]) == C[1]` (reduction of order) when the reduced 1st-order ODE
can't be closed. 0 FAIL (A(R)=0 symbolic + numeric gate).

## Done
- [x] μ(x,y′) search (`dsolve_ifactor.c`): Cases A (Kamke 226), C (`(1+y'²)/(x−y)`),
      D (Kamke 66) + Lemma-2 μ̃; `if_factor_select` (FactorList), `if_mutilde`.
- [x] Per-candidate gate `ifactor_gate_R` (build R + A(R)=0 symbolic `ifactor_R_ok`
      + numeric `ifactor_R_num_ok`); shared `ifactor_find_R`; explicit path unchanged.
- [x] Emit: `dsolve_run_first_integral` + `dsolve_method_builtin_first_integral`
      (dsolve_common); decline IVP; d/dx(R)|y''→Φ verify. Cascade AFTER lie2.
- [x] Pinned builtin `DSolve`ReducibleFirstIntegral` (ATTR_PROTECTED + docstring).
- [x] Tests: `t_m56_*` (test_dsolve.c, verified in isolation) + stress
      `test_dsolve_m56_stress.c` (A/C/D grids) — PASSES.
- [x] Stress suites m5/m12/m14/m18/m55/m56 green; Stage-1 (Airy) preserved;
      check-c99 green; loop smoke-test stable.
- [x] Docs: calculus.md (2 entries), changelog 2026-09-21.md, DSOLVE_PLAN.md M56
      block + M18 Stage 2 resolved. version.h 0.171→0.172.
- [x] Corpus run 1 (under load): 585 P / 619 U / **0 FAIL**; +17 vs baseline 568.
      24 U→P (~20 `_reducible _mu_*` = M56; confirmed 1156 μ(x,y'), 198/1157 μ(x,y)
      first integrals, all sympy=False); 7 P→U = documented `_with_linear_symmetries`
      timing-flaky cluster (M56 declines linear ODEs → not regressions).

## Remaining (blocked on clean corpus re-run number)
- [ ] Regenerate reports/2.1.2.{tsv,md} from clean run
- [ ] STATUS.md M56 scoreboard row + bucket table
- [ ] Fill corpus number in changelog + DSOLVE_PLAN.md; lower CMakeLists gate baseline
- [ ] Commit + tag v0.172

## Pre-existing (NOT M56 — flagged to user)
- dsolve_tests aborts at t_m19_whittaker_confluent (SpecialFunctionForm), in-suite only
  (passes in isolation). Proven pre-existing by clean-tree stash rebuild. Blocks the
  whole unit-suite tail (M20–M56) via ASSERT exit(1). Out of scope; worth its own task.

## Review
M56 complete (v0.172). μ(x,y′) integrating factors (Lemma-3 A/C/D + Lemma-2 μ̃) emitted as
reduction-of-order first integrals `R(x,y[x],y'[x])==C[1]` via the new `dsolve_run_first_integral`
runner + `DSolve`ReducibleFirstIntegral`. Correctness = A(R)=0 symbolic + numeric gate, so a
wrong μ only declines. §2.1.2 **568 → 585 (+17), 0 FAIL**, all in the 2nd_reducible_mu bucket
(38→55); the 7 P→U are the pre-existing 8 s timing cluster (M56 declines linear). Explicit
Stage-1 path unchanged (Airy preserved). Stress m5/m12/m14/m18/m55/m56 green; check-c99 green.
Deferred: Stage 3 μ(y,y′), general Case D, Cases E/F. Pre-existing dsolve_tests t_m19 in-suite
abort flagged (not M56).
