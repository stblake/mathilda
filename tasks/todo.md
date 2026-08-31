# DSolve continuation — 1a FirstOrderSubstitution + 1d AutonomousReduction

Continuing `DSOLVE_PLAN.md`. Both new methods produce solutions verifiable by
back-substitution (so the substrate verify + tests genuinely check them) and fill
gaps the current cascade declines (confirmed by REPL probe 2026-08-31).

## Phase 1 — DSolve`FirstOrderSubstitution (1a)
- [x] `src/calculus/dsolve_fos.c`: `y'[x] == F(a x + b y + c)`.
  - F = solve_top_derivative(P,1); FY = F with y[x]->Y symbol.
  - r = D[FY,x]/D[FY,Y]; require r free of {x, Y}, D[FY,Y] != 0, D[FY,x] != 0.
  - H(W) = FY /. Y -> (W - r x); require free of x (genuine combination).
  - Solve autonomous v' == r + H(v) inline: ∫dW/(r+H) == x + C[1], Solve for W.
  - y = Wsol - r x, one branch per Solve root.
- [x] Register in dsolve.c cascade (after reduce_order) + method enum/string.
- [x] Unit tests: (x+y)^2, (x+y+1)^2, Sqrt[x+y], (2y-x)^2, named-method, IVP.
- [x] Stress test: randomized F(a x + b y + c) family, back-sub verified.

## Phase 2 — DSolve`AutonomousReduction (1d)
- [x] `src/calculus/dsolve_autonomous.c`: `y''[x] == f(y, y')` missing x.
  - F = solve_top_derivative(P,2); Ftmp = F with y'[x]->P, y[x]->Y symbols.
  - require ds_free_of(Ftmp, x) (autonomous) AND contains Y (distinguishes from
    ReductionOfOrder's missing-y case).
  - Stage 1: p p'(Y) == Ftmp[P->p[Y]] ; recurse DSolve for p(Y); freeze C[1].
  - Stage 2: y' == p(y) ; recurse DSolve; rename frozen const -> C[2].
  - Guard: final body must depend on x (reject degenerate y=const).
- [x] Register in dsolve.c cascade (after fos) + method enum/string.
- [x] Unit tests: y y''==(y')^2, y''==(y')^2/y, named-method.
- [x] Stress test: randomized autonomous family, back-sub verified.

## Phase 3 — gates + docs
- [x] `tests/test_dsolve.c`: add unit + stress cases; ctest green.
- [x] valgrind delta check (baseline-subtract per lessons); make check-c99.
- [x] docs/spec/builtins/calculus.md: add both rows (+ ReductionOfOrder row).
- [x] docs/spec/changelog/2026-08-31.md: new weekly file + summary.
- [x] Mathilda_spec.md changelog table row.
- [x] REPL spot-checks: DSolve[...], ?DSolve`FirstOrderSubstitution, etc.

## Review

Both methods landed, all gates green.

- **DSolve`FirstOrderSubstitution** (`src/calculus/dsolve_fos.c`) — `y'==F(a x + b
  y + c)` via constant ratio `r = F_x/F_y`, `v = y + r x`, inline autonomous
  separation. Solves `(x+y)^2`, `(x+y+1)^2`, `(2y-x)^2`, `(m x+y)^2`. Declines
  cleanly when the antiderivative does not invert (`Sqrt[x+y]`, `(x+2y)^3`).
- **DSolve`AutonomousReduction** (`src/calculus/dsolve_autonomous.c`) — `y''==f(y,
  y')` missing x, two-stage recursion (`p p'(y)==f`, then `y'==p(y)`) with
  cross-stage constant renumbering and an x-dependence guard against the
  degenerate `y=const`. Solves `a y y'' == b (y')^2` family → exp/power/reciprocal
  forms. Declines on elliptic (`y''==2y^3`, `y''==-Sin[y]`).
- **Cascade** — appended after `ReductionOfOrder`; both added to the `Method ->`
  enum/string and `DSolve`<Method>[...]` registrars.
- **Tests** — 11 new unit cases + 2 randomized stress families (24 solve/verify
  pairs) in `tests/test_dsolve.c`; `check_solves` asserts `Head === List` before
  the residual (a declined DSolve gives a vacuous `PossibleZeroQ True`). All
  `dsolve_tests` pass.
- **Gates** — `make check-c99` clean; valgrind shows no leak delta over a
  baseline running the same sub-engines (13.4 KB vs 14.3 KB definitely-lost — the
  pre-existing shutdown noise).
- **Docs** — `docs/spec/builtins/calculus.md` (table + examples, incl. the
  previously-missing `ReductionOfOrder` row), `docs/spec/changelog/2026-08-31.md`
  (new weekly file), `Mathilda_spec.md` changelog row, `DSOLVE_PLAN.md` 1a/1d
  items marked done, `tasks/lessons.md` (HoldAll-variable + vacuous-verify traps).

### Notes / follow-ups
- `DSolve` HoldAll **removed** (follow-up request): attributes are now
  `Protected` (matching Mathematica), applied to DSolve and every
  `DSolve`<Method>`. Equations stored in a variable now solve. All 62 dsolve_tests
  pass unchanged; `t_not_holdall` locks it in.
- FirstOrderSubstitution's IVP constant-fit is weak for transcendental (Tan/Log)
  solutions where `Solve` cannot invert the condition — left as general solution;
  not a regression.
- EnergyIntegral is now subsumed by AutonomousReduction for elementary cases;
  a dedicated elliptic/`WeierstrassP` path remains future work.
