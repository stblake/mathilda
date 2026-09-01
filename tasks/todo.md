# DSolve: fix noisy failures on Chini / integrating-factor first-order ODEs

## Tasks
- [x] 1. Guard ordering: reduce `nexp` (Cancel) before free-of guard (dsolve_bernoulli.c, dsolve_chini.c)
- [x] 2. Chini second reduction (linear-term removal → separable) in dsolve_chini_first_integral
- [x] 3. Exact via x^a y^b integrating factor in dsolve_exact.c (gated to PolynomialQ M,N)
- [x] 4. Message hygiene: mute speculative arithmetic across builtin_dsolve
- [x] 5. Build + verify both examples solve cleanly (no Power::infy)
- [x] 6. Regression: pure Bernoulli, autonomous Chini, bare exact, elliptic decline stay clean
- [x] 7. DSolve unit + stress tests + check-c99 + valgrind
- [x] 8. Docs: spec/builtins DSolve entry, changelog 2026-08-31, lessons.md, memory

## Review

**Outcome.** Both reported ODEs now solve; DSolve no longer leaks
`Power::infy`/`Infinity::indet`.
- Ex1 `y'==-x E^-x - y + x E^(2x) y^3` → implicit `∫du/(u^3-1) - x^2/2 == C[1]`,
  `u = y[x] E^x` (Chini reduction b). Verified: numeric residual ~4e-16.
- Ex2 `(x y-2x)y'==y-y^2+3x^2 y^3` → two explicit branches of
  `3x + 1/(xy) - 1/(xy^2) == C[1]` (exact via `μ=x^-2 y^-3`). Verified: ~4e-16.

**Root causes (three defects).**
1. Bernoulli/Chini tested the exponent `nexp` for FreeQ before reducing it; the
   unexpanded `Q` sent `ds_free_of` to the numeric sampler, which hit `0^(-1)` and
   misread the resulting ComplexInfinity as "non-zero" → valid equation rejected +
   leaked message. Fixed by reducing with **Cancel** before the guard.
2. Chini only did the reduce-to-autonomous reduction (a). Added reduction (b),
   linear-term removal → separable.
3. Exact only searched `μ(x)`/`μ(y)`. Added `μ=x^a y^b`, gated to polynomial M,N.
Plus: muted speculative arithmetic warnings across the whole DSolve cascade.

**Regression caught and fixed during work.** First attempt used `Simplify` (not
Cancel) to reduce `nexp`, and ran the `x^a y^b` `Solve` ungated. Both blew up on
the radical systems from `AutonomousReduction` (`y''==2y^3` → `y'==Sqrt[y^4+K]`),
turning a 27.9s suite into a >60s `alarm()` timeout. Fixed with `Cancel` +
`PolynomialQ` gate; suite back to 28.4s, all pass.

**Gates.** `dsolve_tests` / `dsolve_stress_tests` / `dsolve_m5_stress_tests` all
green (2 new cases: `t_chini_linremoval`, `t_exact_xayb`); `make check-c99` clean;
valgrind — no new definitely/possibly-lost over the ~13.6KB one-time engine
baseline.

**Files.** `dsolve_bernoulli.c`, `dsolve_chini.c`, `dsolve_exact.c`, `dsolve.c`,
`tests/test_dsolve.c`; docs `docs/spec/builtins/calculus.md` +
`docs/spec/changelog/2026-08-31.md`.
