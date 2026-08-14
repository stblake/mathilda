# NMinimize RandomSearch: compiled local polish

## Root cause (verified)
`NMinimize[..., Method -> {"RandomSearch", "SearchPoints" -> 1000, "Method" -> "NelderMead"}]`
took 24.8s (≈10–25× Mathematica). Isolated:
- DifferentialEvolution 0.98s, NelderMead(global) 0.59s, RandomSearch(1000) 24.1s.
- RandomSearch scales linearly with SearchPoints (~25 ms/start).
- Profile dominated by the interpreter: `fm_eval_with_bindings`, `apply_own_values`,
  `add_rule`, `symtab_add_own_value`, `builtin_plus`, `builtin_abs`, `expr_copy/free`.

The compiled objective (`NmDriver.f_prog`) is consumed **only** by global-search
scoring (`nm_eval_obj`). The **local polish** (`nm_local_polish` → `fm_run_bfgs`,
findmin.c:3540) passes `D->f_raw` and evaluates through `fm_eval_scalar` /
`fm_eval_with_bindings` — 20 OwnValue installs + deep copy + full evaluate +
numericalize per point. RandomSearch does one polish per SearchPoint ⇒ 1000
interpreter BFGS runs. (The objective itself already compiles; `Compile[]` already
lowers `Sum`. The gap is the polish, not the objective.)

## Fix — compiled objective + FD gradient in the local solvers (broad)
Route every local-solver objective evaluation through the compiled program, and
take the gradient by finite differences off it. Elegant, minimal surface, general
(all NMinimize methods), interpreter fallback preserved, answers unchanged.

- [ ] Add file statics `g_fm_obj_expr` / `g_fm_obj_prog` / `g_fm_obj_nargs`.
- [ ] Fast path at top of `fm_eval_scalar`: when `f == g_fm_obj_expr` and arity
      matches and `compiled_eval_real` returns finite, use it; else interpreter.
      Pointer identity keeps constraint/gradient sub-exprs on the interpreter.
- [ ] `nm_minimize_driver`: set the statics (save/restore) around the solve,
      pointing at `f_eff` / `f_prog`; restore before `compiled_free(f_prog)`.
- [ ] Force FD gradient when compiled: pass `NULL` for `g_exprs` in the polish
      calls (`nm_local_polish`, `nm_continuous_solve`) when `D->f_prog` set, so
      `fm_grad_finite_diff` (→ `fm_eval_scalar` → compiled) is used instead of the
      interpreter symbolic-gradient path.

## Verify
- [ ] Rebuild; the reported example drops from ~25s toward ~1s, same optimum.
- [ ] `tests/` full suite (esp. test_nminimize, test_findmin) — no regressions.
- [ ] MPFR final refinement still exact (uses symbolic gradient, untouched).
- [ ] Docs + changelog (week of 2026-08-10).

## Notes / follow-ups
- FindMinimum uses a separate driver (`findmin_driver`) with no `f_prog`; the
  static is inert there → unchanged. Compiling there is a clean follow-up.
- General (non-box) constraints in the local penalty polish still evaluate on the
  interpreter (different Expr pointers); follow-up via `g_progs`.

## Review (done 2026-08-14)
Root cause was NOT the objective failing to compile — it already compiles, and
`Compile[]` already lowers `Sum`. The compiled program (`f_prog`) served only the
global-search *scoring*; the local polish (`fm_run_bfgs`) ran on the interpreter,
and RandomSearch does one polish per SearchPoint. Fix routes the polish's
objective evals through `f_prog` via a pointer-identity fast path in
`fm_eval_scalar`, with an FD gradient taken off the compiled objective.

Implemented (all in `src/numerical_calculus/findmin.c`):
- statics `g_fm_obj_expr/_prog/_nargs` + fast path in `fm_eval_scalar`.
- register/deregister in `nm_minimize_driver` (plain reset, not save/restore — a
  compilable objective cannot contain a nested NMinimize, so no active
  registration is ever clobbered; also dodges a GCC-16 `-Wmaybe-uninitialized`).
- `nm_local_polish` / `nm_continuous_solve` pass `g_exprs = NULL` when `f_prog`
  set → FD-on-compiled gradient.

Results: RandomSearch(1000) 24.8 s → 0.54 s (identical `1.76e-5`); DE 0.98 → 0.023;
NelderMead 0.59 → 0.13. Clean `-Wall -Wextra` build, `make check-c99` clean,
`nminimize_tests` + `findmin_tests` all pass. Docs + changelog updated.
