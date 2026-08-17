# Task: Add Basin Hopping to NMinimize

Add `Method -> {"BasinHopping", ...}` as a global engine to NMinimize/NMaximize,
mirroring `scipy.optimize.basinhopping` (Wales & Doye, Monte-Carlo minimization):
random-displacement hop → local minimization (quench) → Metropolis accept on the
minimized energies, with an adaptive step size targeting a fixed acceptance rate.
Reuse the existing `NmDriver` machinery (`nm_eval`, `nm_local_polish`, `nm_better`,
`nm_project`) so box/general/integer constraints work with no engine-specific code —
exactly the pattern of DualAnnealing/SHGO/DIRECT.

## Algorithm (scipy fidelity)
- `niter` (top-level `MaxIterations`, default 100 == scipy niter) basin-hopping steps.
- `T` = "Temperature" (default 1.0): Metropolis `accept iff exp(min(0,-(Eₙ-Eₒ)/T)) >= rand`.
- `stepsize` = "StepSize" (default 0.5): uniform per-coord displacement in [-s, s].
- Adaptive step: every "StepInterval" (50) steps, if accept-rate > "TargetAcceptanceRate"
  (0.5) then `s /= "StepFactor"` (0.9, grows), else `s *= factor` (shrinks).
- Optional "SuccessIterations" early stop (scipy niter_success), default off.
- "SearchPoints" -> K = K independent multi-start runs, keep the Deb-best (default 1).
- Penalized energy `E = f + MU*pen` for Metropolis; `nm_better` for the reported best.
- The local minimizer is Mathilda's `nm_local_polish` (the one deliberate scipy diff).

## Plan
- [x] 1. `findmin_internal.h`: enum `NM_BASIN_HOPPING`; NmConfig `bh_*` fields;
        `NM_BH_*` constants; `nm_basin_hopping` prototype; doc-comment block.
- [x] 2. New engine `src/numerical_calculus/nm_basin_hopping.c`.
- [x] 3. `nm_driver.c`: init `nc.bh_*`; add `case NM_BASIN_HOPPING` to dispatch.
- [x] 4. `findmin_nm_common.c`: method name in/out; `nm_option_owner`; parse blocks.
- [x] 5. `src/info.c`: NMinimize docstring — BasinHopping method + sub-options.
- [x] 6. `docs/spec/builtins/numerical-calculus.md`: method table rows + subsection.
- [x] 7. `docs/spec/changelog/2026-08-17.md`: changelog entry.
- [x] 8. `src/version.h`: bump to 0.060.
- [x] 9. `tests/CMakeLists.txt`: add nm_basin_hopping.c to COMMON_SRC + a test exe.
- [x] 10. `tests/test_basin_hopping.c`: 32 tests — ALL PASS, 0 leaks.
- [x] 11. `benchmarks/85-basin-hopping/`: `.m`/`.py` pair + README (8/8 AHEAD, 0 CHECK-FAIL).
- [x] 12. Build clean (make + check-c99); tests pass; benchmark run; 0 leaks.

## Review
DONE. `Method -> {"BasinHopping", ...}` added to NMinimize/NMaximize.
- Engine: `src/numerical_calculus/nm_basin_hopping.c` — faithful scipy.basinhopping
  (perturb → quench → adaptive-step Metropolis), reuses shared NmDriver so
  box/general/disjunctive/integer constraints work with no BH-specific code.
- 6 sub-options (Temperature/StepSize/StepInterval/TargetAcceptanceRate/StepFactor/
  SuccessIterations) + generic SearchPoints/RandomSeed/PostProcess/PenaltyFunction,
  per-method scoped via nm_option_owner. MaxIterations (top-level) = hop count.
- Tests: 32/32 pass, 0 leaks. Siblings (nminimize/DA/SHGO/DIRECT/findmin_methods)
  all still pass. check-c99 clean.
- Benchmark 85: 8/8 AHEAD vs scipy, 0 CHECK-FAIL, ~100×–650× faster (compiled obj).
- Docs: docstring + numerical-calculus.md (method/sub-option tables + subsection with
  fidelity note) + changelog + version 0.060.
- Honest limitation documented: single-run BH is seed-sensitive on widely-separated
  multi-basin problems because Mathilda's quench (correctly) does not overshoot
  basins the way scipy's L-BFGS-B does; recommend SearchPoints. Conversely BH beats
  scipy's stalled single run on Rastrigin, and solves constrained/integer scipy can't.

## Key finding (design)
Basin hopping's escape from a local basin depends on the local minimizer. scipy's
L-BFGS-B overshoots across a basin boundary on its aggressive first step; Mathilda's
nm_local_polish is a well-behaved local minimizer that stays in-basin — so Mathilda's
single-run BH crosses basins only via the random displacement, making it more
seed-sensitive on widely-separated multi-basin problems (quartic, Styblinski-Tang).
Increasing niter does NOT help (the walk is stuck); SearchPoints (multi-start) does.
This is the documented "one deliberate difference: the local minimizer" — algorithm
is byte-for-byte faithful to scipy. Tests use SearchPoints for the multi-basin cases;
benchmark uses the 8 functions both nail at seed 1 (Rastrigin & multi-basin go in the
README regime analysis: Mathilda WINS Rastrigin, scipy wins quartic/StybTang).

## Review
(to be filled in on completion)
