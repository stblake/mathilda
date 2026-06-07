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

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
