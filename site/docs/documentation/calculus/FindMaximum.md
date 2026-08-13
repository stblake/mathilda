# FindMaximum

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindMaximum[f, {x, x0}]`**

searches for a local maximum of f starting from x = x0.

**`FindMaximum[f, {x, x0, x1}]`**

derivative-free 1D search bracketing the maximum from two starts (Brent on -f).

**`FindMaximum[f, {x, xstart, xmin, xmax}]`**

bracketed 1D Brent search on \[xmin, xmax\] starting from xstart.

**`FindMaximum[f, {{x, x0}, {y, y0}, ...}]`**

n-D local maximum from a user-supplied start.

**`FindMaximum[f, {x, y, ...}]`**

n-D local maximum auto-starting each variable at 0.

**`FindMaximum[{f, cons}, vars]`**

local maximum subject to box and Inequality constraints.

<details>
<summary>Notes</summary>

Methods (Method -\> ...): Automatic           picks Brent for 1D, QuasiNewton (BFGS) for n-D. "Brent"             derivative-free golden-section + parabolic interpolation; 1D only; honours MPFR WorkingPrecision. "QuasiNewton"       BFGS with cubic line search; uses analytic Gradient if given, otherwise central differences; default for n\>=2; honours MPFR WorkingPrecision. "ConjugateGradient" Polak-Ribiere CG with line search; lower memory than BFGS for large n; gradient-based. "Newton"            full Hessian step via modified Cholesky factorization; falls back to a steepest-descent step when the Hessian is not negative definite or unavailable. Options: Method              algorithm selector (see above). WorkingPrecision    MachinePrecision (double) or a positive digit count (MPFR; honoured by Brent and BFGS). MaxIterations       positive integer cap on outer iterations; default 500. AccuracyGoal        Automatic | Infinity | digits; absolute tolerance on |f| (and |x| where applicable). PrecisionGoal       Automatic | Infinity | digits; relative tolerance on step size. Gradient            Automatic (finite differences) or an explicit list { dfdx1, dfdx2, ... } in the same order as vars.  The gradient is taken with respect to f, not -f. StepMonitor         :\> body run after each accepted step, with the variables locally bound to their current values. EvaluationMonitor   :\> body run on every function/gradient evaluation. FindMaximum has HoldAll and effectively uses Block to localize the variables.  Internally maximises by minimising -f, then negates the objective value in the result.  Returns {fmax, {x -\> xmax, ...}}.

</details>

## Examples (9)

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

### Applications (3)

```mathematica
In[7]:= FindMaximum[Sin[x], {x, 1}]
Out[7]= {1.0, {x -> 1.5708}}

In[8]:= FindMaximum[x (10 - x), {x, 0}]
Out[8]= {25.0, {x -> 5.0}}

In[9]:= FindMaximum[Sin[x] Sin[2 y], {{x, 1}, {y, 1}}]
Out[9]= {1.0, {x -> 1.5708, y -> 0.785398}}
```

## Algorithm

findmin.c

FindMinimum / FindMaximum — Mathematica-compatible local numerical optimization. Both have HoldAll | Protected attributes and use a Block-style snapshot/restore of the search variables' OwnValues so that user-level definitions of those names are not perturbed during iteration.

Supported forms ---------------

```text
  FindMinimum[f,           {x, x0}]                  1D, Brent default
  FindMinimum[f,           {x, x0, x1}]              1D, two-start bracket
  FindMinimum[f,           {x, xstart, xmin, xmax}]  1D, bracket
  FindMinimum[f,           {{x, x0}, {y, y0}, ...}]  n-D, QuasiNewton default
  FindMinimum[f,           {x, y, ...}]              n-D, auto start = 1
  FindMinimum[{f, cons},   vars]                     constrained
```

Options (Rule[...] in trailing position, any order):

```text
  Method            -> Automatic | "Brent" | "Newton" | "QuasiNewton"
                                 | "ConjugateGradient"
  WorkingPrecision  -> MachinePrecision | digits   (MPFR for Brent + BFGS)
  MaxIterations     -> positive integer (default 500)
  AccuracyGoal      -> Automatic | Infinity | digits
  PrecisionGoal     -> Automatic | Infinity | digits
  Gradient          -> Automatic | { dfdx1, dfdx2, ... }
  StepMonitor       -> :> body
  EvaluationMonitor -> :> body
```

Constraints (inside the {f, cons} form): boolean tree of comparisons.

```text
  Box  ( a <= x <= b , x >= a , x <= b , etc. on a bare variable )
    → enforced by projection after each iterate.
  General ( g(x) <= 0 , h(x) == 0 , etc. )
    → quadratic-penalty wrapper around the inner solver; outer μ schedule.
  Or[...] / Element / Integers → emit FindMinimum::nimpl and return NULL.
```

Output: { f_min, { x -> x_min, y -> y_min, ... } }. FindMaximum returns { f_max, ... } via a thin wrapper that minimises −f and negates the first component of the result.

Returns NULL (unevaluated) on any failure — variable bindings are always restored to their pre-call state, even on the error path.

## Implementation notes

**Algorithm.** `FindMaximum` (`HoldAll | Protected`) is a thin wrapper over
`FindMinimum` (src/findmin.c, `builtin_findmaximum`): it negates the objective,
runs the same local optimizer, and negates the first component of the resulting
`{f_min, {x -> x_min, ...}}` pair to report `{f_max, {x -> x_max, ...}}`. All
machinery — Brent in 1-D, BFGS quasi-Newton / conjugate-gradient / Newton in
n-D, symbolic gradients/Hessian with a central-difference fallback, Armijo line
search, box-projection and quadratic-penalty constraint handling, MPFR extended
precision — is inherited unchanged from `FindMinimum`. See `FindMinimum` for the
full description.

**Complexity / limits.** Same as `FindMinimum`: local search only. The negation
is precision-aware (`mpfr_neg` for `EXPR_MPFR` results, plain real otherwise).

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [FindMinimum](../../calculus/FindMinimum/), [Block](../../scoping-constructs/Block/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [PrecisionGoal](../../other-advanced/PrecisionGoal/)

- J. Nocedal, S. J. Wright, *Numerical Optimization*, 2nd ed. (Springer, 2006).
- Source: [`src/findmin.c`](https://github.com/stblake/mathilda/blob/main/src/findmin.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_findmin.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin.c)

## Notes & additional examples

### Notes

`FindMaximum[f, {x, x0}]` returns `{fmax, {x -> xmax, ...}}`. Internally it
maximises by minimising `-f`, so the same Brent (1D) and BFGS quasi-Newton
(n-D) machinery as `FindMinimum` applies. The first example recovers the
peak of `Sin` at `x = π/2`; the multivariate case locates a saddle-free
maximum of the product `Sin[x] Sin[2 y]` at `(π/2, π/4)`.
