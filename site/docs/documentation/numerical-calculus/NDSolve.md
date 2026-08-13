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

## Examples (7)

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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [InterpolatingFunction](../../functional-programming/InterpolatingFunction/), [Derivative](../../calculus/Derivative/), [Dt](../../calculus/Dt/), [Interpolation](../../functional-programming/Interpolation/), [PrecisionGoal](../../other-advanced/PrecisionGoal/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [ComplexExpand](../../arithmetic/ComplexExpand/)

- Source: [`src/numerical_calculus/ndsolve.c`](https://github.com/stblake/mathilda/blob/main/src/numerical_calculus/ndsolve.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_ndsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndsolve.c)
- Tests: [`tests/test_ndsolve_classical.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndsolve_classical.c)
- Tests: [`tests/test_ndsolve_pde.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndsolve_pde.c)
