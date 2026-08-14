# NMinimize: DifferentialEvolution options + 20-D rotated Rastrigin stress test

## Requests
1. Add the 20-D rotated Rastrigin as an NMinimize unit test.
2. Ensure all `"DifferentialEvolution"` sub-options work as expected.
3. Use Storn & Price's 10n = 200 default population.

## Findings (verified before coding)
- All DE options were *parsed*, but two were consumed only by NelderMead:
  - `"Tolerance"` → used only in `nm_neldermead` (simplex spread); DE's early-break
    used the AccuracyGoal/PrecisionGoal-derived threshold. **Silently ignored by DE.**
  - `"InitialPoints"` → used only in `nm_simplex_from_points`; DE always
    random-initialised. **Silently ignored by DE.**
  - `"CrossProbability"`/`"ScalingFactor"` accepted any value with no validation.
- Explicit `"DifferentialEvolution"` clamped the population to `[15, 40]`, so a
  20-D solve ran 40 members instead of 10n = 200.

## Changes (all in src/numerical_calculus/findmin.c)
- [x] Population = `Clip[10·d, {15, 200}]` for explicit DE and Automatic alike.
- [x] DE convergence break honors an explicit `"Tolerance"` (default path unchanged).
- [x] New `nm_population_from_points` seeds the leading DE population members from
      `"InitialPoints"`, rest random; RNG stream unchanged when absent.
- [x] Validate `"CrossProbability"` ∈ [0,1] and `"ScalingFactor"` ∈ (0,2]; invalid
      warns `NMinimize::sopt` and falls back (matches the SA sub-options).

## Tests (tests/test_nminimize.c)
- [x] `test_rotated_rastrigin_stress` — loose shape/feasibility/basin invariant.
- [x] `test_de_options_effective` — Tolerance, InitialPoints, CR, F each steer the
      search; invalid CR/F warn + fall back.

## Docs
- [x] docs/spec/builtins/numerical-calculus.md (population, Tolerance, InitialPoints,
      CR/F validation).
- [x] docs/spec/changelog/2026-08-10.md (this ISO week).

## Review / results
- Full `nminimize_tests` (63 tests) pass; `findmin_tests` pass; `make check-c99` clean.
- No regressions: `test_griewank_differentialevolution` (n=10, now NP=100) still
  passes; existing seeded runs with no InitialPoints reproduce (RNG stream intact).
- Stress test lands at 37.8084 (feasible, 20 coords) — a good basin, not the global
  0. Plain DE does not crack 20-D rotated Rastrigin; the bound is a "found a good
  region" check, robust to cross-platform QR + DE drift.
- Empirically confirmed each option now takes effect: Tolerance 5.5e-16 (tight) vs
  3.5e-3 (loose); InitialPoints-at-origin → raw 0.0 vs random 65.
