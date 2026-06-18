# NSolve stress-test suite

## Goal
Create and run a comprehensive stress-test suite for `NSolve`, exercising every
option and both machine and arbitrary precision, beyond the functional coverage
already in `tests/test_nsolve.c`.

## Probed behaviour (ground truth, 2026-06-18)
- `x^n` polynomials: full root count, Reals filter correct, residuals tiny.
- Wilkinson literal `(x-1)...(x-10)`: 10 roots, accurate to <1e-2, imag ~0.
- Complex coefficients (`x^4+(2+I)x+(1-3I)`): full count, tiny residuals.
- Arbitrary precision: residual scales with WorkingPrecision (30/50/80);
  `Precision[root]` ~= requested digits.
- `MaxRoots` sweep caps correctly (0,1,5,12 -> exact; >count -> count).
- `Method` "EndomorphismMatrix"/"Homotopy"/"Symbolic" + Automatic all agree.
- `RandomSeeding`: seed CHANGES how many of the N solutions the randomized
  eigenvalue engine recovers (seed 1 -> 5, seed 999 -> 6 on a 6-solution system).
  => seed tests assert VALIDITY of returned solutions, not a fixed count.
- `VerifySolutions` True/False/Automatic all return valid solutions.
- `PrecisionGoal`, `MaxIterations` accepted (no-op on poly engine).
- Domains: Reals filters, Integers via Solve, Complexes default.
- Transcendental (`Cos[x]==x`, `E^x-x==7`) via FindRoot grid seeding.
- Known edge (NOT asserted): `NSolve[{eqn},{x}]` list-wrapped single univariate
  equation routes to Solve which itself bails -> stays unevaluated.

## Plan
- [x] Probe real NSolve behaviour to ground expectations
- [x] Write `tests/test_nsolve_stress.c`
- [x] Register `nsolve_stress_tests` in `tests/CMakeLists.txt`
- [x] Build and run; all pass
- [x] valgrind (diff vs baseline noise)
- [x] Update changelog

## Review
- `tests/test_nsolve_stress.c` adds 17 stress test functions; all pass.
  Existing `nsolve_tests` still passes (test-only change, no NSolve edits).
- Leak check done on the MAIN binary (the test binary's `alarm(60)` from
  test_utils.h fires under valgrind's ~30x slowdown and kills the run mid-eval,
  reporting in-flight allocations as "lost" — not a real signal). NSolve
  workload adds 0 definitely-lost bytes over the 12,800 B / 400 block macOS
  dyld baseline, across poly / system / transcendental / Integers paths.
- Notable characterizations captured as tests:
  - `RandomSeeding` changes the recovered solution SUBSET of a system (seed 1
    -> 5 of 6, seed 999 -> 6) => assert validity, not count, for seeds.
  - `Product[x-k,{k,..}]` does NOT expand here -> Wilkinson uses literal factors.
  - Edge NOT asserted: `NSolve[{eqn},{x}]` (list-wrapped single univariate eqn)
    routes to Solve, which bails -> stays unevaluated. Pre-existing Solve limit.
