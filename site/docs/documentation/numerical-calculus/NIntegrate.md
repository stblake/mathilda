# NIntegrate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NIntegrate[f, {x, xmin, xmax}]
    gives a numerical approximation to the integral of f with respect to x from xmin to xmax.

NIntegrate[f, {x, xmin, xmax}, {y, ymin, ymax}, ...] evaluates a multidimensional integral by adaptive cubature over a constant box, or iterated 1D quadrature when an inner bound depends on an outer variable. The variable is localised (HoldAll). xmin/xmax may be Infinity, -Infinity, or complex (a straight-line contour); extra nodes {x, x0, x1, ..., xk} give a piecewise-linear contour or mark interior singularities. Method -> Automatic chooses globally-adaptive Gauss-Kronrod for smooth finite integrands, double-exponential (tanh-sinh / sinh-sinh / exp-sinh) for endpoint singularities and infinite ranges and high precision, a Levin/zeros scheme for oscillatory integrands, an exponential endpoint map plus integration-between-the-zeros for an oscillatory endpoint singularity, and Monte-Carlo for high dimensions and region (Boole) integrands. Machine or arbitrary precision via WorkingPrecision.

Options: Method (Automatic | GlobalAdaptive | GaussKronrodRule | DoubleExponential | TrapezoidalRule | LevinRule | OscillatorySingularity | MonteCarlo | QuasiMonteCarlo | AdaptiveMonteCarlo | PrincipalValue), WorkingPrecision (default MachinePrecision), PrecisionGoal, AccuracyGoal, MaxRecursion, MinRecursion, MaxPoints, Exclusions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NIntegrate[Cos[x], {x, 0, Pi/2}]
Out[1]= 1.0

In[2]:= NIntegrate[Exp[-x^2], {x, 0, Infinity}]
Out[2]= 0.886227

In[3]:= NIntegrate[1/Sqrt[x], {x, 0, 1}]
Out[3]= 2.0

In[4]:= NIntegrate[Sin[x]/x, {x, 0, Infinity}]
Out[4]= 1.5708

In[5]:= NIntegrate[Exp[-x^2 - y^2], {x, -Infinity, Infinity}, {y, -Infinity, Infinity}]
Out[5]= 3.14159
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
