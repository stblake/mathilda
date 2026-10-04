---
source: src/numerical_calculus/ndsolve.c
references:
  - "J. R. Dormand and P. J. Prince, *A family of embedded Runge–Kutta formulae*, J. Comput. Appl. Math. **6** (1980) 19–26."
  - "E. Hairer, S. P. Nørsett and G. Wanner, *Solving Ordinary Differential Equations I: Nonstiff Problems*, 2nd ed. (Springer, 1993) — the initial-step heuristic (§II.4)."
  - "B. Fornberg, *Generation of finite difference formulas on arbitrarily spaced grids*, Math. Comp. **51** (1988) 699–706 — the method-of-lines spatial stencils."
---
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
