### Worked examples

```mathematica
In[1]:= NDSolve[{y'[x] == y[x], y[0] == 1}, y, {x, 0, 1}]  (* the solution is returned as an InterpolatingFunction *)
Out[1]= {{y -> InterpolatingFunction[{{0.0, 1.0}}, <>]}}
```

```mathematica
In[1]:= First[y[1] /. NDSolve[{y'[x] == y[x], y[0] == 1}, y, {x, 0, 1}]]  (* y[1] = E *)
Out[1]= 2.71828
```

```mathematica
In[1]:= First[y[1.5] /. NDSolve[{y''[x] + y[x] == 0, y[0] == 0, y'[0] == 1}, y, {x, 0, 10}]]  (* a second-order ODE: y = Sin[x] *)
Out[1]= 0.997495
```

```mathematica
In[1]:= First[x[1] /. NDSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0}, {x, y}, {t, 0, 10}]]  (* a first-order system: x = Cos[t] *)
Out[1]= 0.540302
```

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
