# Task: Disjunctive (Or) constraints in NMinimize / NMaximize

## Goal
Support constraints like `(x-3)^2+(y-2)^2<=0.1 || (x+2.8)^2+(y+3.1)^2<=0.1`
in NMinimize's global optimizers (RandomSearch/DE/SA/NelderMead), instead of
emitting `NMinimize::nimpl: disjunctive (Or) constraints are not yet supported`.

## Design
- A disjunction is feasible iff at least one branch is feasible ⇒ its penalty is
  the MIN over branches of that branch's violation penalty (0 when any branch
  holds). This composes with the existing Deb feasibility gate (`nm_better`)
  with zero tuning.
- Global (derivative-free) search consumes the min-penalty; the smooth local
  polish (BFGS) ignores disjunctions, and the post-polish `nm_eval` +
  `nm_better` gate rejects any polish that leaves the disjunctive-feasible set.
- FindMinimum (gradient penalty method) still rejects `Or` — a non-smooth min
  breaks its μ-schedule / gradient path. NMinimize-only feature.

## Steps
- [ ] Add `FmDisjunction { Expr* expr; }` struct near `FmGenCon`.
- [ ] Add `fm_bool_supported(Expr*)` structural validator (And/Or/Inequality/cmp).
- [ ] Add `fm_bool_penalty(...)` recursive penalty: And=Σ, Or=min, Inequality=Σ
      pairs, comparison=squared violation (honours "PenaltyFunction").
- [ ] Extend `fm_collect_constraints` with a disjunction sink (NULL ⇒ Or nimpl);
      collect `Or[...]` subtrees when the sink is provided.
- [ ] Thread NULL sink through the FindMinimum call site (Or stays nimpl there).
- [ ] Add `disj`/`ndisj` to `NmDriver`; add disjunction sum to `nm_eval_pen`.
- [ ] Wire collection + free of `disj` into the NMinimize setup/cleanup.
- [ ] Docstring/spec/changelog updates (builtin modified).
- [ ] Unit test: Himmelblau constrained to two disks -> min at (3,2), f≈0.

## Verification
- [ ] Build main + nminimize_tests clean (`-std=c99 -Wall -Wextra`).
- [ ] Run the user's In[5] example: no nimpl warning, returns feasible ≈0 result.
- [ ] `make check-c99` clean.
- [ ] Full nminimize_tests suite passes; new test deterministic.

## Review — DONE (2026-08-14)

Implemented in `src/numerical_calculus/findmin.c`:
- `FmDisjunction` struct + `fm_bool_supported` (structural validator) +
  `fm_bool_penalty` (recursive And=Σ / Or=min / Inequality=Σ / cmp=squared).
- `fm_collect_constraints` gained a disjunction sink; FindMinimum passes NULL
  (Or still `nimpl`), NMinimize passes the real sink.
- `nm_eval_pen` adds `Σ_disj min-branch penalty`.
- Disjunction-aware **local polish** (`nm_polish_gens` + `nm_collect_branch_gens`
  + `fm_free_gens`): folds each disjunction's active (min-penalty) branch into
  the smooth BFGS penalty solve so RandomSearch (pure multi-start polish) doesn't
  strand feasible starts in the infeasible gap.

**Non-obvious bug found & fixed:** symbolic `D[...]` taken *during* the search
returns 0 because the optimisation variables are transiently value-bound
(differentiating a constant). Fixed by using finite-difference gradients for the
polish constraints (`grad_exprs = NULL`); the objective gradient is unaffected
(differentiated once at setup where the variables are free).

Verified:
- User In[5] (Himmelblau two-disks) → `~0` at (3,2), feasible. No `nimpl`.
- `x^2, x<=-2||x>=2` → 4 at x=±2 across DE/SA/NelderMead/**RandomSearch**, with
  and without an enclosing box.
- Boundary optima (nearest-point-to-disk, NMaximize disjunction) correct.
- FindMinimum still rejects `Or` (`nimpl`).
- `findmin_tests` + `nminimize_tests` (incl. new `test_disjunctive_constraints`)
  all pass; `make check-c99` clean; valgrind shows no leak in the new code
  (only the known macOS objc/dyld baseline).

Docs: docstring (info.c), spec (numerical-calculus.md), changelog (2026-08-10.md).
