# NDSolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NDSolve[eqns, u, {x, xmin, xmax}]
    solves the ordinary differential equations eqns numerically for the
    function u over xmin <= x <= xmax, returning {{u -> InterpolatingFunction[...]}}.
NDSolve[eqns, {u1, u2, ...}, {x, xmin, xmax}] solves a system.
NDSolve[eqns, u[x], {x, xmin, xmax}] gives u[x] -> InterpolatingFunction[...][x].
    Higher-order equations (u''[x] == ...) are reduced to first order.
NDSolve[eqns, u, {t, tmin, tmax}, {x, xmin, xmax}] solves a partial
    differential equation over a rectangular region by the method of lines,
    giving a 2-D InterpolatingFunction applied as u[t, x].
    Options: Method, WorkingPrecision, AccuracyGoal, PrecisionGoal,
    MaxSteps, MaxStepSize, MaxStepFraction, StartingStepSize,
    InterpolationOrder, StepMonitor, EvaluationMonitor.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/numerical_calculus/ndsolve.c`](https://github.com/stblake/mathilda/blob/main/src/numerical_calculus/ndsolve.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
