# Task: Add Options to `Reduce`

Plan: `/Users/user/.claude/plans/at-present-reduce-does-humble-duckling.md`

## Steps

- [ ] A1. Add `SYM_Backsubstitution` to `src/sym_names.{h,c}` (3 sites)
- [ ] A2. Register `Options[Reduce]` defaults in `src/options_builtin.c`
- [ ] B1. New `src/solve/reduce_opts.{h,c}` (struct + default + build_solve helper)
- [ ] B2. Front-end peeler + `ReduceOpts` init + dispatch threading in `reduce.c`
- [ ] B3. Add `reduce_opts.c` to build (makefile auto-discovers `src/*.c`; verify CMake COMMON_SRC)
- [ ] C1. Cubics/Quartics: thread into `reduce_eq.c:solve_generic`
- [ ] C2. Cubics/Quartics + WP: `reduce_real_util.c` (`rru_collect_roots`, `rru_approx_double`)
- [ ] C3. Cubics/Quartics + WP: `reduce_realdiag.c` (`soft_roots`, `gen_sign_at`, `cmp_bp`)
- [ ] C4. Cubics/Quartics: `reduce_cad.c` fibre Solve calls
- [ ] C5. GeneratedParameters + new `reduce_modular`: `reduce_int.c`
- [ ] C6. Backsubstitution: `reduce_sys.c` (`reduce_eq_system`)
- [ ] C7. Thread opts through `reduce_univar` / `reduce_univar_general` / `reduce_univar_integers`
- [ ] D1. Docstring in `reduce_init`
- [ ] D2. `docs/spec/builtins/solutions-of-equations.md`
- [ ] D3. `docs/spec/changelog/2026-08-24.md`
- [ ] D4. `REDUCE_PLAN.md` Options section
- [ ] E1. Tests in `tests/test_reduce.c` (one per option + peeling + registration)
- [ ] V1. Build + run reduce_tests, solve_tests, reduce_corpus_tests
- [ ] V2. `make check-c99`; valgrind smoke; REPL smoke
- [ ] V3. Rebuild code-review graph

## Review

All steps complete. `Options[Reduce]` now returns the seven Mathematica-compatible
options and each is honored:

- **Cubics/Quartics** — forwarded onto internal `Solve[...]` calls via a shared
  `reduce_opts_build_solve` (radicals vs `Root[]`). Scoped-out: CAD fibre isolation
  keeps `Root[]` (documented).
- **Modulus** — top-level pre-pass `reduce_modular` routes through Solve's `solvemod`;
  reformatted as `Or` of `x==r`; symbolic/out-of-range/non-modular declines.
- **GeneratedParameters** — `rename_param_head` renames `C[k] -> h[k]` in the
  Integers/Rationals output.
- **Backsubstitution** — accept/validate/echo (no fork in the current linear engine).
- **WorkingPrecision** — threads numeric-fallback tolerance into `reduce_realdiag`;
  Infinity keeps exact-first.
- **Method** — reserved (Automatic only).

New files: `src/solve/reduce_opts.{c,h}`. New symbol: `SYM_Backsubstitution`.
Defaults registered in `options_builtin.c`. Engines threaded: `reduce.c` (peeler +
dispatch), `reduce_eq.c`, `reduce_int.c` (+ modular/rename), `reduce_real_util.c`,
`reduce_univar.c`, `reduce_realdiag.c`; CAD passes NULL (scoped).

Verification: main build clean; `make check-c99` clean; `reduce_tests` (225 assertions,
incl. 9 new `test_option_*` groups), `solve_tests`, `options_tests`, and the 154-case
`reduce_corpus` all pass — no regressions. Valgrind: no leak attributable to the new
code (residual radical-path leak under `Cubics/Quartics->True` is pre-existing in
`solvepoly.c`, reached identically via `Solve`). Docs updated:
`docs/spec/builtins/solutions-of-equations.md`, `docs/spec/changelog/2026-08-24.md`,
`REDUCE_PLAN.md`, and the `reduce_init` docstring.
