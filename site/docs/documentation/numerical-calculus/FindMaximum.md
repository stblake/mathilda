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

Methods (Method -\> ...): Automatic           picks Brent for 1D, QuasiNewton (BFGS) for n-D. "Brent"             derivative-free golden-section + parabolic interpolation; 1D only; honours MPFR WorkingPrecision. "QuasiNewton"       BFGS with cubic line search; uses analytic Gradient if given, otherwise central differences; default for n\>=2; honours MPFR WorkingPrecision. "ConjugateGradient" Polak-Ribiere CG with line search; lower memory than BFGS for large n; gradient-based. "Newton"            full Hessian step via modified Cholesky factorization; falls back to a steepest-descent step when the Hessian is not negative definite or unavailable. "LBFGSB"            limited-memory BFGS (history 10) with box constraints via active-set projection; O(m\*n) per step, so it scales to large n where QuasiNewton's dense Hessian is O(n^2); aliases "LBFGS", "LimitedMemoryBFGS" (a Mathilda extension). "Powell"            derivative-free conjugate directions with Brent line search; no gradient, so it suits non-smooth or black-box objectives; supports box bounds (not general constraints); alias "PrincipalAxis". "NelderMead"        derivative-free downhill simplex; no gradient; robust on smooth black-box objectives (weak on non-smooth -- prefer "Powell" there); supports box bounds (not general constraints). "TNC"               Hessian-free truncated Newton: inner conjugate gradient solving the Newton system via Hessian-vector products (finite differences of the gradient), with active-set box bounds; O(n) memory, so it scales where "Newton" (full symbolic Hessian) cannot, using truer curvature than "LBFGSB"; general constraints via the penalty wrapper; alias "TruncatedNewton". "SLSQP"             sequential least-squares QP (Han-Powell SQP) for smooth constrained problems: each step solves a QP (BFGS Lagrangian-Hessian, linearized constraints) accepted by an L1-penalty line search; handles equality, inequality AND box constraints DIRECTLY (not via the penalty wrapper), so it reaches the true constrained optimum with super-linear local convergence; reduces to damped-BFGS when unconstrained; alias "SequentialQuadraticProgramming". "COBYLA"            Powell's derivative-free linear-approximation trust-region method: models f and every constraint by a linear approximation and solves a two-stage trust-region LP (feasibility then objective) with an L-infinity penalty merit; the only derivative-free method that accepts general (non-box) constraints, so it suits non-smooth / black-box CONSTRAINED objectives; equalities handled by splitting; machine precision only. "COBYQA"            derivative-free trust-region SQP with QUADRATIC interpolation models (vs COBYLA's linear): captures curvature, so it converges tighter on smooth problems and navigates curved valleys (Rosenbrock) that COBYLA cannot; handles equality + inequality + bound constraints natively; machine precision only. "NewtonCG"          line-search truncated Newton: inner CG solves the Newton system via Hessian-vector products with negative-curvature truncation, then a Wolfe line search; Hessian-free, scales to large n; UNCONSTRAINED (rejects constraints, ignores bounds); alias "Newton-CG". "Dogleg"            trust-region Powell dogleg: blends the Cauchy and Newton points on the dense Hessian along the dogleg path; falls back to steepest-to-boundary when the model is not positive definite; UNCONSTRAINED; alias "dogleg". "TrustNCG"          trust-region Steihaug-Toint truncated CG using Hessian-vector products only; Hessian-free; UNCONSTRAINED; aliases "trust-ncg", "TrustRegionNewtonCG". "TrustExact"        trust-region with a near-exact More-Sorensen subproblem on the dense Hessian (Levenberg shift + hard-case branch); handles indefinite Hessians dogleg cannot; UNCONSTRAINED; alias "trust-exact". "TrustKrylov"       trust-region GLTR: Lanczos-tridiagonalize the Hessian in the Krylov space of the gradient (Hessian-vector products), solve the tridiagonal subproblem exactly; Hessian-free, handles the indefinite/hard case; UNCONSTRAINED; alias "trust-krylov". Options: Method              algorithm selector (see above). WorkingPrecision    MachinePrecision (double) or a positive digit count (MPFR; honoured by Brent and BFGS). MaxIterations       positive integer cap on outer iterations; default 500. AccuracyGoal        Automatic | Infinity | digits; absolute tolerance on |f| (and |x| where applicable). PrecisionGoal       Automatic | Infinity | digits; relative tolerance on step size. Gradient            Automatic (finite differences) or an explicit list { dfdx1, dfdx2, ... } in the same order as vars.  The gradient is taken with respect to f, not -f. StepMonitor         :\> body run after each accepted step, with the variables locally bound to their current values. EvaluationMonitor   :\> body run on every function/gradient evaluation. FindMaximum has HoldAll and effectively uses Block to localize the variables.  Internally maximises by minimising -f, then negates the objective value in the result.  Returns {fmax, {x -\> xmax, ...}}.

</details>

## Examples (11)

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

### Applications (4)

```mathematica
In[8]:= FindMaximum[Sin[x], {x, 1}]
Out[8]= {1.0, {x -> 1.5708}}

In[9]:= FindMaximum[x (10 - x), {x, 0}]
Out[9]= {25.0, {x -> 5.0}}

In[10]:= FindMaximum[Sin[x] Sin[2 y], {{x, 1}, {y, 1}}]
Out[10]= {1.0, {x -> 1.5708, y -> 0.785398}}

In[11]:= FindMaximum[10 - (x - 3)^2 - (y + 1)^2, {{x, 0}, {y, 0}}, Method -> "LBFGSB"]
Out[11]= {10.0, {x -> 3.0, y -> -1.0}}
```

## Implementation notes

**Algorithm.** `FindMaximum` (`HoldAll | Protected`) is a thin wrapper over
`FindMinimum` (src/numerical_calculus/findmin.c, `builtin_findmaximum`): it negates the objective,
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

**See also:** [FindMinimum](../../numerical-calculus/FindMinimum/), [Block](../../scoping-constructs/Block/), [NMinimize](../../numerical-calculus/NMinimize/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [PrecisionGoal](../../other-advanced/PrecisionGoal/), [Compile](../../control-flow/Compile/)

- J. Nocedal, S. J. Wright, *Numerical Optimization*, 2nd ed. (Springer, 2006).
- Source: [`src/numerical_calculus/findmin.c`](https://github.com/stblake/mathilda/blob/main/src/numerical_calculus/findmin.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_findmin.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin.c)
- Tests: [`tests/test_findmin_cobyla.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin_cobyla.c)
- Tests: [`tests/test_findmin_cobyqa.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin_cobyqa.c)
- Tests: [`tests/test_findmin_dogleg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findmin_dogleg.c)

## Notes & additional examples

### Notes

`FindMaximum[f, {x, x0}]` returns `{fmax, {x -> xmax, ...}}`. Internally it
maximises by minimising `-f`, so the same Brent (1D) and BFGS quasi-Newton
(n-D) machinery as `FindMinimum` applies. The first example recovers the
peak of `Sin` at `x = π/2`; the multivariate case locates a saddle-free
maximum of the product `Sin[x] Sin[2 y]` at `(π/2, π/4)`.

Every `FindMinimum` method is available, including `Method -> "LBFGSB"`
(limited-memory BFGS with bound constraints; aliases `"LBFGS"`,
`"LimitedMemoryBFGS"`), shown in the last example locating the peak of a
concave paraboloid at `(3, -1)`.
