# FindMinimum

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindMinimum[f, {x, x0}]`**

searches for a local minimum of f starting from x = x0.

**`FindMinimum[f, {x, x0, x1}]`**

derivative-free 1D search bracketing the minimum from two starts (Brent).

**`FindMinimum[f, {x, xstart, xmin, xmax}]`**

bracketed 1D Brent search on \[xmin, xmax\] starting from xstart.

**`FindMinimum[f, {{x, x0}, {y, y0}, ...}]`**

n-D local minimum from a user-supplied start.

**`FindMinimum[f, {x, y, ...}]`**

n-D local minimum auto-starting each variable at 0.

**`FindMinimum[{f, cons}, vars]`**

local minimum subject to box and Inequality constraints.

<details>
<summary>Notes</summary>

Methods (Method -\> ...): Automatic           picks Brent for 1D, QuasiNewton (BFGS) for n-D. "Brent"             derivative-free golden-section + parabolic interpolation; 1D only; honours MPFR WorkingPrecision. "QuasiNewton"       BFGS with cubic line search; uses analytic Gradient if given, otherwise central differences; default for n\>=2; honours MPFR WorkingPrecision. "ConjugateGradient" Polak-Ribiere CG with line search; lower memory than BFGS for large n; gradient-based. "Newton"            full Hessian step via modified Cholesky factorization; falls back to a steepest-descent step when the Hessian is not positive definite or unavailable. "LBFGSB"            limited-memory BFGS (history 10) with box constraints via active-set projection; O(m\*n) per step, so it scales to large n where QuasiNewton's dense Hessian is O(n^2); aliases "LBFGS", "LimitedMemoryBFGS" (a Mathilda extension). "Powell"            derivative-free conjugate directions with Brent line search; no gradient, so it suits non-smooth or black-box objectives; supports box bounds (not general constraints); alias "PrincipalAxis". "NelderMead"        derivative-free downhill simplex; no gradient; robust on smooth black-box objectives (weak on non-smooth -- prefer "Powell" there); supports box bounds (not general constraints). "TNC"               Hessian-free truncated Newton: inner conjugate gradient solving the Newton system via Hessian-vector products (finite differences of the gradient), with active-set box bounds; O(n) memory, so it scales where "Newton" (full symbolic Hessian) cannot, using truer curvature than "LBFGSB"; general constraints via the penalty wrapper; alias "TruncatedNewton". "SLSQP"             sequential least-squares QP (Han-Powell SQP) for smooth constrained problems: each step solves a QP (BFGS Lagrangian-Hessian, linearized constraints) accepted by an L1-penalty line search; handles equality, inequality AND box constraints DIRECTLY (not via the penalty wrapper), so it reaches the true constrained optimum with super-linear local convergence; reduces to damped-BFGS when unconstrained; alias "SequentialQuadraticProgramming". "COBYLA"            Powell's derivative-free linear-approximation trust-region method: models f and every constraint by a linear approximation and solves a two-stage trust-region LP (feasibility then objective) with an L-infinity penalty merit; the only derivative-free method that accepts general (non-box) constraints, so it suits non-smooth / black-box CONSTRAINED objectives; equalities handled by splitting; machine precision only. "COBYQA"            derivative-free trust-region SQP with QUADRATIC interpolation models (vs COBYLA's linear): captures curvature, so it converges tighter on smooth problems and navigates curved valleys (Rosenbrock) that COBYLA cannot; handles equality + inequality + bound constraints natively; machine precision only. "NewtonCG"          line-search truncated Newton: inner CG solves the Newton system via Hessian-vector products with negative-curvature truncation, then a Wolfe line search; Hessian-free, scales to large n; UNCONSTRAINED (rejects constraints, ignores bounds); alias "Newton-CG". "Dogleg"            trust-region Powell dogleg: blends the Cauchy and Newton points on the dense Hessian along the dogleg path; falls back to steepest-to-boundary when the model is not positive definite; UNCONSTRAINED; alias "dogleg". "TrustNCG"          trust-region Steihaug-Toint truncated CG using Hessian-vector products only; Hessian-free; UNCONSTRAINED; aliases "trust-ncg", "TrustRegionNewtonCG". "TrustExact"        trust-region with a near-exact More-Sorensen subproblem on the dense Hessian (Levenberg shift + hard-case branch); handles indefinite Hessians dogleg cannot; UNCONSTRAINED; alias "trust-exact". "TrustKrylov"       trust-region GLTR: Lanczos-tridiagonalize the Hessian in the Krylov space of the gradient (Hessian-vector products), solve the tridiagonal subproblem exactly; Hessian-free, handles the indefinite/hard case; UNCONSTRAINED; alias "trust-krylov". Options: Method              algorithm selector (see above). WorkingPrecision    MachinePrecision (double) or a positive digit count (MPFR; honoured by Brent and BFGS). MaxIterations       positive integer cap on outer iterations; default 500. AccuracyGoal        Automatic | Infinity | digits; absolute tolerance on |f| (and |x| where applicable). PrecisionGoal       Automatic | Infinity | digits; relative tolerance on step size. Gradient            Automatic (finite differences) or an explicit list { dfdx1, dfdx2, ... } in the same order as vars. StepMonitor         :\> body run after each accepted step, with the variables locally bound to their current values. EvaluationMonitor   :\> body run on every function/gradient evaluation. FindMinimum has HoldAll and effectively uses Block to localize the variables.  Returns {fmin, {x -\> xmin, ...}}.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= FindMinimum[(x - 3)^2 + 1, {x, 0}]
Out[1]= {1.0, {x -> 3.0}}

In[2]:= FindMinimum[x Cos[x], {x, 2}]
Out[2]= {-3.28837, {x -> 3.42562}}

In[3]:= FindMinimum[x Cos[x], {x, 7, 1, 15}]
Out[3]= {-9.47729, {x -> 9.52933}}

In[4]:= FindMinimum[Sin[x] Sin[2 y], {{x, 2}, {y, 2}}]
Out[4]= {-1.0, {x -> 1.5708, y -> 2.35619}}

In[5]:= FindMinimum[(1-x)^2 + 100 (y-x^2)^2, {{x, 0}, {y, 0}}]
Out[5]= {1.58322e-20, {x -> 1.0, y -> 1.0}}

In[6]:= FindMaximum[Cos[x], {x, 0}]
Out[6]= {1.0, {x -> -2.3206e-09}}
```

### Options (1)

```mathematica
In[7]:= FindMinimum[(x - Pi)^2, {x, 0}, WorkingPrecision -> 50]
Out[7]= {0.0, {x -> 3.1415926535897932384626433832795028841971693993751}}
```

### Applications (6)

```mathematica
In[8]:= FindMinimum[x^2 - 4 x + 7, {x, 0}]
Out[8]= {3.0, {x -> 2.0}}

In[9]:= FindMinimum[Cos[x] + x/5, {x, 0, 10}]
Out[9]= {-0.391749, {x -> 2.94023}}

In[10]:= FindMinimum[(1 - x)^2 + 100 (y - x^2)^2, {{x, -1}, {y, 1}}]
Out[10]= {3.46541e-23, {x -> 1.0, y -> 1.0}}

In[11]:= FindMinimum[Gamma[x], {x, 1.5}]
Out[11]= {0.885603, {x -> 1.46163}}

In[12]:= FindMinimum[(1 - x)^2 + 100 (y - x^2)^2, {{x, -1.2}, {y, 1}}, Method -> "LBFGSB"]
Out[12]= {4.71192e-22, {x -> 1.0, y -> 1.0}}

In[13]:= FindMinimum[(x - 2)^2 + (y - 3)^2, {{x, 0, 0, 1}, {y, 0, 0, 1}}, Method -> "LBFGSB"]
Out[13]= {5.0, {x -> 1.0, y -> 1.0}}
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| FindMinimum Rosenbrock 2-D | 0.19 s | 0.158 s | 3.57 s |
| Fit degree-5 polynomial, 500 points | 0.096 s | 0.076 s | 0.038 s |
| FindMinimum quadratic | 0.016 s | 0.061 s | 0.172 s |
| Fit linear, 200 points | 0.01 s | 0.015 s | 0.015 s |
| FindMinimum x Sin[x] on [2,6] | 0.005 s | 0.297 s | 1.14 s |

## Implementation notes

**Algorithm.** `FindMinimum` (`HoldAll | Protected`) performs local numerical
optimization (src/numerical_calculus/findmin.c), Block-snapshotting and restoring the search
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

## References

**See also:** [FindMaximum](../../numerical-calculus/FindMaximum/), [Block](../../scoping-constructs/Block/), [NMinimize](../../numerical-calculus/NMinimize/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [PrecisionGoal](../../other-advanced/PrecisionGoal/), [Compile](../../control-flow/Compile/)

- J. Nocedal, S. J. Wright, *Numerical Optimization*, 2nd ed. (Springer, 2006) — BFGS, line search, conjugate gradient.
- W. H. Press et al., *Numerical Recipes*, 3rd ed. (Cambridge, 2007) — Brent's method.
- Source: [`src/numerical_calculus/findmin.c`](https://github.com/stblake/mathilda/blob/main/src/numerical_calculus/findmin.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_findmin.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin.c)
- Tests: [`tests/test_findmin_cobyla.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin_cobyla.c)
- Tests: [`tests/test_findmin_cobyqa.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin_cobyqa.c)
- Tests: [`tests/test_findmin_dogleg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin_dogleg.c)

## Notes & additional examples

### Notes

`FindMinimum[f, {x, x0}]` performs a local search from the start `x0`,
returning `{fmin, {x -> xmin, ...}}`. The third example is the notorious
Rosenbrock banana valley: BFGS quasi-Newton drives the iterate into the
curved trough and locates the global minimum `(1, 1)` to machine precision.
The Gamma example finds the minimum of the Gamma function on the positive
axis (a root of the digamma function) at `x ≈ 1.4616`.

`Method -> "LBFGSB"` selects **limited-memory BFGS with bound constraints**
(aliases `"LBFGS"`, `"LimitedMemoryBFGS"`; a Mathilda extension). It keeps only
the last 10 correction pairs — `O(m·n)` per step rather than the full-memory
QuasiNewton's `O(n²)` — so it scales to large variable counts, and it handles
box bounds by an active-set projection. The two `LBFGSB` examples above show an
unconstrained solve (the Rosenbrock valley, reached to machine precision) and a
box-constrained solve whose unconstrained minimum `(2, 3)` lies outside the unit
box, so the constrained optimum is pinned to the corner `(1, 1)` with value `5`.
