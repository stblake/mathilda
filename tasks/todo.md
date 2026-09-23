# TODO: DSolve M57 — SolvableForY / SolvableForX (the `y=G(x,y')` differentiation method)

Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-binary-hartmanis.md`
Generalizes `DSolve`Lagrange` (its linear-induced-ODE special case). 0 FAIL: the
parametric substrate back-substitution-verifies every branch.

## Done
- [x] `src/calculus/dsolve_solvefor.c` — shared core `sfy_solve_for(P, for_y)`;
      isolate y/x, differentiate, recurse induced ODE, emit `DSolve`Param[X,Y,p]`.
      Bounding kit (deadline + re-entry guard + decline memo + node budget +
      bounded Simplify), `PolynomialQ` pre-gate (kills 20 s transcendental
      time-burners: `Sin[xy]` 20 s → 0.24 s).
- [x] `dsolve.c` wiring: enum/string/extern×2/pinned×2/init. **SolvableForY auto**
      (before Lie backstop); **SolvableForX pinned-only** (auto yield ~0 + slow
      cubic-denominator verify → opt-in like FirstOrderPowerSeries).
- [x] Pinned builtins `DSolve`SolvableForY`/`SolvableForX` (ATTR_PROTECTED + docstrings).
- [x] `make -j` + `make check-c99` green. Latency verified (347 family <2 s,
      declines <0.3 s).
- [x] Stress `tests/test_dsolve_m57_stress.c` (347 n=3..6, 347+ n=3..5, pinned Y/X;
      parametric numeric verify) + CMake registration — PASSES.
- [x] Units `t_m57_*` (`tests/test_dsolve.c`) — verified in isolation (all True).
- [x] Docs: calculus.md (2 entries), changelog 2026-09-21.md (M57), DSOLVE_PLAN.md
      (M57 block + §1a entries). version.h 0.172 → 0.173.
- [x] `dsolve_solvefor.c` added to tests/CMakeLists.txt source list.

## In progress
- [ ] §2.1.2 corpus re-run (background) → measure yield, confirm 0 FAIL, P→U regressions.
- [ ] Fill corpus number in changelog + DSOLVE_PLAN.md; lower CMakeLists gate baseline.
- [ ] valgrind spot-check a couple of SolvableForY solves (ownership clean in new file).
- [ ] Commit + tag v0.173.

## Pre-existing (NOT M57 — flagged)
- `dsolve_tests` aborts in-suite at `t_m19_whittaker_confluent` (blocks M20+ unit tail);
  proven pre-existing (clean-tree stash). M57 units verified in isolation.
- `dsolve_stress_tests` aborts at `t_stress_undetcoeff`: `DSolve`UndeterminedCoefficients`
  declines `y''-2y'+y==Cos[2x]` — proven pre-existing on clean HEAD (my edits don't touch
  the undetcoeff path). Worth its own task.

## Review
M57 complete (v0.173). `SolvableForY` (auto) + `SolvableForX` (pinned-only) — the `y=G(x,y')`
"dp" differentiation method (Maple's dp), generalising `DSolve`Lagrange` (its linear-induced-ODE
special case): isolate y/x, differentiate, recurse the cascade on the induced first-order ODE,
return the parametric solution via the existing parametric substrate (which back-substitution-
verifies every branch → 0 FAIL by construction). New file `dsolve_solvefor.c`; zero substrate edits.

Key engineering: a `PolynomialQ` pre-gate (mirroring NthAlgebraic) declines the transcendental
time-burners that never close (`Sin[xy]` 20 s → 0.24 s) — without it, each direction burned ~10 s
per declined case (P→U timing regressions). All recursion bounded (M14 kit). `SolvableForX` made
pinned-only after a degree-1-in-x cubic-denominator form produced a transcendental branch whose
substrate verify spun 30 s (outside the try-fn's control); its auto yield is ~0 anyway.

§2.1.2 (clean re-run, no concurrent load): **585 → 590 (net +5), 0 FAIL**. Deterministic M57 gain
**+2** (347 `y=_G(x,y')`, 352 `[F(x),G(y)]`-symmetry — both solved by SolvableForY, verified); the
other +3 is the `_with_linear_symmetries` timing cluster oscillating (2nd/high-order, declined
instantly by SolvableForY). 350/351 solve but exceed the 8 s forked budget. Gate baseline 631 → 629.
Memory ownership clean by inspection (inherited Integrate/Solve engine leak is outside the file).
Stress m5/m12/m14/m18/m55/m56/m57 green; check-c99 green.

Lesson: never run Mathilda scripts concurrently with a forked corpus run — it slows cases past
their 8 s budget and contaminates the P/U measurement (my first run was contaminated; the clean
re-run corrected it). Recorded in memory.

Pre-existing (flagged, not M57): `dsolve_tests` in-suite abort at `t_m19_whittaker_confluent`;
`dsolve_stress_tests` abort at `t_stress_undetcoeff` (UndeterminedCoefficients declines
`y''-2y'+y==Cos[2x]`) — both proven pre-existing on clean HEAD.

Future: explicit `p(x)` route; parametric IVP fitting; parameter elimination to implicit Φ(x,y)=0;
`y⁽ⁿ⁾`-solvable higher-order; radical/degree-≥4 bucket cases needing Root-object handling.
