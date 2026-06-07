# FindMinimum

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FindMinimum[f, {x, x0}]
    searches for a local minimum of f starting from x = x0.
FindMinimum[f, {x, x0, x1}]
    derivative-free 1D search bracketing the minimum from two starts (Brent).
FindMinimum[f, {x, xstart, xmin, xmax}]
    bracketed 1D Brent search on [xmin, xmax] starting from xstart.
FindMinimum[f, {{x, x0}, {y, y0}, ...}]
    n-D local minimum from a user-supplied start.
FindMinimum[f, {x, y, ...}]
    n-D local minimum auto-starting each variable at 0.
FindMinimum[{f, cons}, vars]
    local minimum subject to box and Inequality constraints.

Methods (Method -> ...):
    Automatic           picks Brent for 1D, QuasiNewton (BFGS) for n-D.
    "Brent"             derivative-free golden-section + parabolic interpolation; 1D only; honours MPFR WorkingPrecision.
    "QuasiNewton"       BFGS with cubic line search; uses analytic Gradient if given, otherwise central differences; default for n>=2; honours MPFR WorkingPrecision.
    "ConjugateGradient" Polak-Ribiere CG with line search; lower memory than BFGS for large n; gradient-based.
    "Newton"            full Hessian step via modified Cholesky factorization; falls back to a steepest-descent step when the Hessian is not positive definite or unavailable.

Options:
    Method              algorithm selector (see above).
    WorkingPrecision    MachinePrecision (double) or a positive digit count (MPFR; honoured by Brent and BFGS).
    MaxIterations       positive integer cap on outer iterations; default 500.
    AccuracyGoal        Automatic | Infinity | digits; absolute tolerance on |f| (and |x| where applicable).
    PrecisionGoal       Automatic | Infinity | digits; relative tolerance on step size.
    Gradient            Automatic (finite differences) or an explicit list { dfdx1, dfdx2, ... } in the same order as vars.
    StepMonitor         :> body run after each accepted step, with the variables locally bound to their current values.
    EvaluationMonitor   :> body run on every function/gradient evaluation.

FindMinimum has HoldAll and effectively uses Block to localize the variables.  Returns {fmin, {x -> xmin, ...}}.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FindMinimum[(x - 3)^2 + 1, {x, 0}]
Out[1]= {1.0, {x -> 3.0}}

In[2]:= FindMinimum[x Cos[x], {x, 2}]
Out[2]= {-3.28837, {x -> 3.42562}}

In[3]:= FindMinimum[x Cos[x], {x, 7, 1, 15}]
Out[3]= {-9.47729, {x -> 9.52933}}

In[4]:= FindMinimum[Sin[x] Sin[2 y], {{x, 2}, {y, 2}}]
Out[4]= {-1.0, {x -> 1.5708, y -> 2.35619}}

In[5]:= FindMinimum[{x Cos[x], 1 <= x && x <= 15}, {x, 7}]
Out[5]= {-9.47729, {x -> 9.52933}}

In[6]:= FindMinimum[(1-x)^2 + 100 (y-x^2)^2, {{x, 0}, {y, 0}}]
Out[6]= {1.58322e-20, {x -> 1.0, y -> 1.0}}

In[7]:= FindMaximum[Cos[x], {x, 0}]
Out[7]= {1.0, {x -> -2.3206e-09}}

In[8]:= FindMinimum[(x - 3)^2, {x, 0}, Method -> "ConjugateGradient"]
Out[8]= {0.0, {x -> 3.0}}
```

## Implementation notes

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

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- J. Nocedal, S. J. Wright, *Numerical Optimization*, 2nd ed. (Springer, 2006) — BFGS, line search, conjugate gradient.
- W. H. Press et al., *Numerical Recipes*, 3rd ed. (Cambridge, 2007) — Brent's method.
- Source: [`src/findmin.c`](https://github.com/stblake/mathilda/blob/main/src/findmin.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
