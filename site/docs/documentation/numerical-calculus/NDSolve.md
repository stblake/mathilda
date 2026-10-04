# NDSolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NDSolve[eqns, u, {x, xmin, xmax}]`**

solves the ordinary differential equations eqns numerically for the function u over xmin \<= x \<= xmax, returning {{u -\> InterpolatingFunction\[...\]}}.

**`NDSolve[eqns, {u1, u2, ...}, {x, xmin, xmax}] solves a system.`**

**`NDSolve[eqns, u[x], {x, xmin, xmax}] gives u[x] -> InterpolatingFunction[...][x].`**

Higher-order equations (u''\[x\] == ...) are reduced to first order.

**`NDSolve[eqns, u, {t, tmin, tmax}, {x, xmin, xmax}] solves a partial`**

differential equation over a rectangular region by the method of lines, giving a 2-D InterpolatingFunction applied as u\[t, x\]. Options: Method, WorkingPrecision, AccuracyGoal, PrecisionGoal, MaxSteps, MaxStepSize, MaxStepFraction, StartingStepSize, InterpolationOrder, StepMonitor, EvaluationMonitor.

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= sol = NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 5}]; y[1] /. sol
Out[1]= {0.367879}
```

= Cos[3]

```mathematica
In[2]:= NDSolve[{y''[x] + y[x] == 0, y[0] == 1, y'[0] == 0}, y, {x, 0, 6}]; y[3.0] /. %
Out[2]= y[3.0]
```

Circle: {Cos t, -Sin t}

```mathematica
In[3]:= NDSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0}, {x, y}, {t, 0, 6}]
Out[3]= {{x -> InterpolatingFunction[{{0.0, 6.0}}, <>], y -> InterpolatingFunction[{{0.0, 6.0}}, <>]}}
```

Wave equation u_tt = u_xx (default adaptive DOPRI5)

```mathematica
In[4]:= NDSolve[{D[u[t, x], {t, 2}] == D[u[t, x], {x, 2}], u[0, x] == Sin[Pi x], Derivative[1, 0][u][0, x] == 0, u[t, 0] == 0, u[t, 1] == 0}, u, {t, 0, 0.5}, {x, 0, 1}]
Out[4]= {{u -> InterpolatingFunction[{{0.0, 0.5}, {0.0, 1.0}}, <>]}}
```

### Options (3)

Stiff

```mathematica
In[5]:= NDSolve[{y'[x] == -1000 (y[x] - Cos[x]) - Sin[x], y[0] == 1}, y, {x, 0, 3}, Method -> "BackwardEuler"]
Out[5]= {{y -> InterpolatingFunction[{{0.0, 3.0}}, <>]}}
```

```mathematica
In[6]:= NDSolve[{y'[x] == y[x], y[0] == 1}, y, {x, 0, 1}, WorkingPrecision -> 30, PrecisionGoal -> 22, MaxSteps -> 200000]
Out[6]= {{y -> InterpolatingFunction[{{0.0, 1.0}}, <>]}}
```

Heat equation u_t = u_xx, Dirichlet, method of lines

```mathematica
In[7]:= sol = NDSolve[{D[u[t, x], t] == D[u[t, x], {x, 2}], u[0, x] == Sin[Pi x], u[t, 0] == 0, u[t, 1] == 0}, u, {t, 0, 0.05}, {x, 0, 1}, Method -> "BDF"]; u[0.05, 0.5] /. sol
Out[7]= {0.610498}
```

### Applications (4)

The solution is returned as an InterpolatingFunction

```mathematica
In[8]:= NDSolve[{y'[x] == y[x], y[0] == 1}, y, {x, 0, 1}]
Out[8]= {{y -> InterpolatingFunction[{{0.0, 1.0}}, <>]}}
```

Y[1] = E

```mathematica
In[9]:= First[y[1] /. NDSolve[{y'[x] == y[x], y[0] == 1}, y, {x, 0, 1}]]
Out[9]= 2.71828
```

A second-order ODE: y = Sin[x]

```mathematica
In[10]:= First[y[1.5] /. NDSolve[{y''[x] + y[x] == 0, y[0] == 0, y'[0] == 1}, y, {x, 0, 10}]]
Out[10]= 0.997495
```

A first-order system: x = Cos[t]

```mathematica
In[11]:= First[x[1] /. NDSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0}, {x, y}, {t, 0, 10}]]
Out[11]= 0.540302
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| NDS Van der Pol mu=1000 | 4.02 s | 0.234 s | 0.362 s |
| NDS harmonic t=100 | 0.884 s | 0.475 s | 37.2 s |
| NDS Van der Pol mu=100 | 0.647 s | 0.271 s | 0.452 s |
| NDSolve Van der Pol mu=10 | 0.441 s | 0.566 s | 3.32 s |
| NDS Lorenz t=5 | 0.337 s | 0.45 s | 14.4 s |
| NDS Van der Pol mu=10 | 0.295 s | 0.559 s | 3.32 s |

## Implementation notes

**Algorithm.** `builtin_ndsolve` → `ndsolve_core` parses `NDSolve[eqns, funcs,
{x, a, b}]`. `nd_scan` walks each equation matching `u[x]`/`Derivative[m][u][x]`
to find every function's maximum order, then **reduces to a first-order system**
`dY/dt = f(t, Y)`: each function of order `n` contributes `n` state components,
chain equations `f[base+m] = y_{base+m+1}`, and a top equation solved
symbolically for the highest derivative (assuming linearity in it, `a·P + b == 0
⇒ P = −b/a`). Non-ODE equations `Derivative[m][u][pt] == v` seed the initial
state `Y0`; a complex RHS/IC is *realified* (Re/Im split, dimension doubled).
Two or more range specs hand off to the method of lines.

The subsystem is factored on two axes: a problem compiled into an `NdProblem`,
and a numerical method exposed as a single-step `NdStepper` driven by one shared
adaptive loop `nd_integrate` (forward to `tmax` and backward to `tmin` from
`t0`). Dispatch is by **name**, not stiffness detection (`nd_lookup_stepper`);
the eight registered steppers:

- **Dormand–Prince 5(4)** (`ndsolve_rk.c`), the Automatic default: a 7-stage
  FSAL embedded explicit Runge–Kutta pair (6 RHS evals per accepted step,
  matching `ode45`), with the embedded `b − b*` estimate driving step control.
- **RK4** (classical fixed-step), **ExplicitEuler**, **ExplicitMidpoint**
  (`ndsolve_euler.c`) — fixed steppers adapted by step-doubling with local
  extrapolation.
- **Adams** (`ndsolve_adams.c`): order-2 PECE (AB2 predictor + trapezoidal
  corrector), RK4 self-start, Milne corrector−predictor error estimate.
- **BackwardEuler**, **ImplicitTrapezoid**, **BDF** (`ndsolve_implicit.c`): the
  stiff path, stage equations solved by Newton iteration (`nd_newton_theta`,
  factoring `I − coef·J` via LAPACK banded/dense LU). BDF is a full
  variable-step variable-order (orders 1–5) engine with exact nonuniform-mesh
  Lagrange coefficients, a Milne predictor estimate, and order ramping.
  `"StiffnessSwitching"` is aliased to BDF.
- **Method of lines** (`ndsolve_mol.c`, `ndsolve_stencil.c`) for a PDE:
  a uniform spatial grid, **Fornberg finite-difference stencils**, Dirichlet/
  Neumann/Robin/periodic BCs, centred or upwind/Lax–Friedrichs advection, the
  big ODE system handed back to the shared driver.

The result is always an `InterpolatingFunction` rule list, built by the generic
`Interpolation` builtin from **cubic-Hermite `{{t_i}, y_i, y'_i}` triples**
(dense output); complex systems recombine paired components as `ifRe + I ifIm`.

**Data structures.** Machine state vectors are flat `double*` of the reduced
dimension `d`; stage scratch is `stages·d`. The RHS has three tiers: a compiled
constant linear operator `A·Y + s(t)` via banded BLAS (`ndsolve_operator.c`,
used as the exact Jacobian when the system is affine); nonlinear bytecode lazily
built through the general `Compile[]` engine (`ndsolve_compile.c`, with a
Curtis–Powell–Reid colored finite-difference Jacobian); and a symbolic sampler
fallback that binds `t` and the state symbols Block-style. Arbitrary precision
(`ndsolve_mpfr.c`) carries `mpfr_t` state/time/step at `out_bits + 64` guard
bits — explicit methods run DOPRI5/RK4 in MPFR, implicit requests an MPFR BDF
(with the scalar step-ratio control done in double). The growable `NdSolution`
holds `ts`, `Ys`, `dYs`.

**Complexity / limits.** Step control is an elementary `0.9·err^(−1/q)`
controller (no PID); the initial step uses Hairer's HNW II.4 heuristic. Options:
`Method`, `WorkingPrecision` (MPFR above machine), `AccuracyGoal` (default
`MachinePrecision`) and `PrecisionGoal` → `rtol`/`atol`, `MaxSteps` (Automatic
backstop 2 000 000; an explicit value is a hard wall), `MaxStepSize`,
`MaxStepFraction` (default 1/10), `StartingStepSize`, `StepMonitor`,
`EvaluationMonitor` (disables the compiled fast paths). `InterpolationOrder` and
`NormFunction` are parsed but unused — interpolation is always cubic Hermite.
NDSolve declines (returns `NULL`) when it cannot identify a well-posed system,
on non-numeric endpoints, or on insufficient initial/boundary conditions;
integration failures (`ndsz`/`ndcf`/`nrnum`/maxsteps) report through
`mth_message` and return the partial solution where one exists. `NDSolve` is
`HoldAll`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [InterpolatingFunction](../../functional-programming/InterpolatingFunction/), [Derivative](../../calculus/Derivative/), [Dt](../../calculus/Dt/), [Interpolation](../../functional-programming/Interpolation/), [PrecisionGoal](../../other-advanced/PrecisionGoal/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [ComplexExpand](../../arithmetic/ComplexExpand/)

- J. R. Dormand and P. J. Prince, *A family of embedded Runge–Kutta formulae*, J. Comput. Appl. Math. **6** (1980) 19–26.
- E. Hairer, S. P. Nørsett and G. Wanner, *Solving Ordinary Differential Equations I: Nonstiff Problems*, 2nd ed. (Springer, 1993) — the initial-step heuristic (§II.4).
- B. Fornberg, *Generation of finite difference formulas on arbitrarily spaced grids*, Math. Comp. **51** (1988) 699–706 — the method-of-lines spatial stencils.
- Source: [`src/numerical_calculus/ndsolve.c`](https://github.com/stblake/mathilda/blob/main/src/numerical_calculus/ndsolve.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_ndsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndsolve.c)
- Tests: [`tests/test_ndsolve_classical.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndsolve_classical.c)
- Tests: [`tests/test_ndsolve_pde.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndsolve_pde.c)

## Notes & additional examples

### Notes

`NDSolve[eqns, funcs, {x, xmin, xmax}]` integrates an ODE initial-value problem
and returns a rule list `{{y -> InterpolatingFunction[...]}}`. The
`InterpolatingFunction` is dense cubic-Hermite output from the `{t, y, y'}` triples
collected along the integration, so it may be evaluated at any point of the
solved interval — **at a machine number**, e.g. `y[1.5]`, not an exact symbolic
argument like `y[Pi]`.

Higher-order ODEs and systems are reduced to a first-order system internally
(each function's chain of derivatives becomes extra state components), so the
second-order and two-function examples above work directly. Complex-valued ODEs
are realified (Re/Im split) and recombined on output.

The default `Method -> Automatic` is **adaptive Dormand-Prince 5(4)** (`"DOPRI5"`,
the embedded explicit Runge-Kutta pair, FSAL, with error-controlled step size);
`"RungeKutta"` is fixed-step RK4, `"ExplicitEuler"`/`"ExplicitMidpoint"` are
adapted by step-doubling, and stiff problems use `Method -> "BDF"` (variable-step
variable-order backward differentiation, orders 1-5) or `"StiffnessSwitching"`,
`"BackwardEuler"`, `"ImplicitTrapezoid"`, `"Adams"`. There is no automatic
stiffness detection — the stiff solver is opt-in. Two or more range specs select
the method of lines for a PDE.

`AccuracyGoal` (default `MachinePrecision`) and `PrecisionGoal` set the local
error tolerance; `WorkingPrecision -> d` runs the explicit/BDF integrators in MPFR
at `d` digits. `MaxSteps`, `MaxStepSize`, and `StartingStepSize` bound the step
loop. `NDSolve` has `HoldAll`.
