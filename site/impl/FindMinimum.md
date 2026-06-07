---
references:
  - "J. Nocedal, S. J. Wright, *Numerical Optimization*, 2nd ed. (Springer, 2006) — BFGS, line search, conjugate gradient."
  - "W. H. Press et al., *Numerical Recipes*, 3rd ed. (Cambridge, 2007) — Brent's method."
source: src/findmin.c
---
**Algorithm.** `FindMinimum` (`HoldAll | Protected`) performs local numerical
optimization (src/findmin.c), Block-snapshotting and restoring the search
variables' OwnValues around the iteration. The variable spec and dimension
choose the inner solver:

- 1-D (`{x, x0}` / two-start / bracket) → **Brent** golden-section/parabolic
  minimization (the default in one dimension).
- n-D (`{{x, x0}, ...}` or bare variable list) → **BFGS quasi-Newton** by
  default, with `"ConjugateGradient"` and a full-Hessian `"Newton"` also
  selectable via `Method`.

Gradients are computed symbolically as a list of `D[f, x_i]`
(`fm_compute_gradient`), with a central-difference numeric fallback when the
symbolic gradient fails (`fm_numeric_gradient`); the Hessian for Newton is the
symbolic `D[D[f, x_i], x_j]` array (`fm_compute_hessian`). The quasi-Newton loop
maintains an approximate inverse Hessian updated by the BFGS formula and takes
steps satisfying an **Armijo backtracking line search**. Constraints in the
`{f, cons}` form are classified: box constraints on bare variables are enforced
by projection after each iterate; general inequality/equality constraints are
handled by a **quadratic-penalty** wrapper around the inner solver with an outer
μ schedule; `Or[...]`/`Element`/`Integers` are rejected with `FindMinimum::nimpl`.

When `WorkingPrecision` requests extended precision and MPFR is built in, the
Brent and BFGS paths run at the requested bit width (`fm_run_bfgs_mpfr`,
`fm_line_search_mpfr`, `fm_eval_gradient_mpfr`); the MPFR BFGS path does not yet
support the penalty/constraint machinery.

`FindMaximum` (`builtin_findmaximum`) is a thin wrapper that minimizes `−f` and
negates the first component of the `{f_min, {x -> ...}}` result.

**Data structures.** `double` arrays for machine-precision gradient/Hessian and
the inverse-Hessian matrix; `mpfr_t` arrays for extended precision. Function and
gradient evaluation re-enter the Mathilda evaluator with current bindings.

**Complexity / limits.** Local minimization only — no global search. BFGS is the
default n-D method (superlinear local convergence); options: `MaxIterations`
(default 500), `AccuracyGoal`/`PrecisionGoal`, user `Gradient`, and held
`StepMonitor`/`EvaluationMonitor`. Returns NULL/unevaluated on non-numeric
evaluation or non-convergence, always restoring variable bindings.
