# FindRoot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FindRoot[f, {x, x0}]
    searches for a numerical root of f starting from x = x0.
FindRoot[lhs == rhs, {x, x0}]
    searches for a numerical solution to the equation.
FindRoot[f, {x, x0, x1}]
    uses a variant of the secant method with x0 and x1 as the first two approximations.
FindRoot[f, {x, xstart, xmin, xmax}]
    uses Brent's method on the bracket [xmin, xmax].
FindRoot[{f1, f2, ...}, {{x, x0}, {y, y0}, ...}]
    searches for a simultaneous numerical root of the system.

Options: Method ('Newton' | 'Secant' | 'Brent' | Automatic), WorkingPrecision, MaxIterations, AccuracyGoal, PrecisionGoal, DampingFactor, Jacobian, StepMonitor, EvaluationMonitor.  FindRoot has HoldAll and effectively uses Block to localize variables.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FindRoot[Sin[x] + Exp[x], {x, 0}]
Out[1]= {x -> -0.588533}

In[2]:= FindRoot[Cos[x] == x, {x, 0}]
Out[2]= {x -> 0.739085}

In[3]:= FindRoot[{y == Exp[x], x + y == 2}, {{x, 1}, {y, 1}}]
Out[3]= {x -> 0.442854, y -> 1.55715}

In[4]:= FindRoot[Sin[x], {x, 3}, WorkingPrecision -> 50]
Out[4]= {x -> 3.14159265358979323846264338328757519744320987120168}

In[5]:= FindRoot[(Cos[z + I] - 2) (z + 2), {z, 1.0 + 0.1 I}]
Out[5]= {z -> -1.66935e-13 + 0.316958*I}

In[6]:= FindRoot[Cos[x] - x, {x, 0, 1}, Method -> "Brent"]
Out[6]= {x -> 0.739085}

In[7]:= FindRoot[x^2 - 2, {x, 1.0, 2.0}, Method -> "Secant"]
Out[7]= {x -> 1.41421}

In[8]:= FindRoot[(x - 1)^3, {x, 0.5}, DampingFactor -> 3]
Out[8]= {x -> 1.0}
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
