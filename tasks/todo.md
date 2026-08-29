# Task: Exponential-tower dilogarithm engine (`rt_cherry_dilog_exp`)

Close `Integrate[x/(-1+E^x), x]` and the exp-tower dilog family via a native
Cherry engine (mirror of `cherry_dilog.c` transported to θ=E^(cx)).

## Plan
`/Users/user/.claude/plans/mathilda-s-implementation-of-cherry-s-temporal-stearns.md`

## Steps
- [x] 1. Prototype the full algorithm as a Mathilda `.m` script; validate on the
      whole corpus (rational×x, outer-log, mixtures, reducible den, declines).
      DONE — `tasks/dilog_exp_prototype.m`. Flagship + family + outer-log +
      mixtures all diff-back 0. `Log[1-E^x]` (negative-arg outer log) declines
      cleanly (reversed-pair/iπ; out of scope, sound).
- [ ] 2. Write `src/calculus/cherry_dilog_exp.{c,h}` transcribing the validated
      prototype (mirror cherry_dilog.c structure + memory discipline).
- [ ] 3. Register in `risch_special.c` (`#include` + `RT_SPECIAL_FORMS[]` row).
- [ ] 4. Build (`make -j`); fix any C99/portability issues (`make check-c99`).
- [ ] 5. Verify the flagship + family via `-file` scripts (diff-back = 0).
- [ ] 6. Write `tests/test_cherry_dilog_exp.c`; wire into `tests/CMakeLists.txt`
      (COMMON_SRC + target); build & run.
- [ ] 7. Regression: cherry_dilog/li/ei + risch/integrate test targets.
- [ ] 8. Docs: builtin page + changelog `docs/spec/changelog/2026-08-25.md`.
- [ ] 9. Update memory / lessons.
- [ ] SIDE TASK (user-requested): benchmark the `.m` prototype vs the C engine
      on the corpus; report the speedup.

## Review

**Done (2026-08-29).** New engine `rt_cherry_dilog_exp` (`src/calculus/cherry_dilog_exp.c`),
the exp-tower mirror of `rt_cherry_dilog`, closes `Integrate[x/(-1+E^x),x] =
x Log[1-E^-x] - PolyLog[2,E^-x]` and the family (`x/(1±E^x)`, `x E^x/(E^x-1)`,
`x/(E^(2x)-1)`, `Log[1±E^x]`, `Log[1+E^-x]`, mixtures). Native Cherry tower matching
(no substitution): ansatz + linear solve over `{θ,x,u_k}` + PowerExpand diff-back.
Registered as a `PolyLog`/`RT_SF_TOP_EXP` form (`risch_special.c`).

- Steps 1–9 complete. Extras done at user request: (a) all five Cherry engines exposed as
  `Integrate`Cherry`{Ei,ExpMultiterm,Li,Dilog,DilogExp}` debug surfaces (hub
  `cherry_builtins_init` in `cherry_driver.c`), each with docstring + PROTECTED;
  (b) speed comparison — C engine ~8% faster than the `.m` prototype for raw matching
  (16.6 vs 18.0 ms/integral), because both share the same C symbolic kernels.
- Key fix: top-level depth gate (`g_integrate_depth > 1 → NULL`) so a DerivativeDivides
  `u=Log[x]` substitution doesn't reroute `Log[x]/(1-x)` through the exp engine and change
  its clean form. See `tasks/lessons.md` (2026-08-29).
- Tests: `cherry_dilog_exp_tests` (new) + all cherry/risch/integrate suites green (39+).
  `make check-c99` clean.
- Docs: `docs/spec/builtins/calculus.md` + `docs/spec/changelog/2026-08-24.md`. Memory:
  `project_cherry_dilog_exp_engine.md`.

**Follow-on (2026-08-29) — general polylogarithm ladder.** `rt_cherry_polylog_exp`
(`src/calculus/cherry_polylog_exp.c`) generalises the exp-tower dilog to ARBITRARY
WEIGHT and ALGEBRAIC ROOTS: `P(x)/Q(E^(cx))` → partial-fraction over roots of Q
(rational/algebraic) + exact Cherry ladder `∫x^n/(θ-ρ) = Σ -(1/ρ)(n!/(n-k)!)/c^{k+1}
x^{n-k} PolyLog[k+1,ρ/θ]`, polylogs up to weight n+1. Closes `x^2/(E^x-1)`,
`x^4/(E^(5x)-1)`, `x/(E^(2x)+E^x-1)`, `x^2/(E^(2x)+E^x-1)` (Q(√5)) — all diff-back 0.
Registered ahead of dilog_exp (cleaner rational forms). Debug surface
`Integrate\`Cherry\`PolyLogExp`. Tests `cherry_polylog_exp_tests`. Updated a stale
dilog_exp decline assertion (x^2/(E^x-1) now integrates). All cherry + risch +
integrate suites (26) green; `make check-c99` clean. NOT yet committed.
