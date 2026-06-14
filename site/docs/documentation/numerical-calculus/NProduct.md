# NProduct

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NProduct[f, {i, imin, imax}]
    gives a numerical approximation to the product of f for i from imin to imax.

NProduct[f, {i, imin, imax, di}] uses step di. imax may be Infinity. NProduct[f, {i, ...}, {j, ...}, ...] evaluates a multidimensional product (an inner bound may depend on an outer index). The index is localised (HoldAll). Evaluated as Exp[NSum[Log[f], ...]], so the NSum engine (Euler-Maclaurin for monotone factors, Wynn's epsilon otherwise) and its convergence test carry over. Machine or arbitrary precision via WorkingPrecision.

Options: Method (Automatic | EulerMaclaurin | WynnEpsilon), WorkingPrecision (default MachinePrecision), NProductFactors (leading factors taken explicitly, default 15), NProductExtraFactors, WynnDegree, VerifyConvergence (default True; a divergent product gives ComplexInfinity), AccuracyGoal, PrecisionGoal.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NProduct[1 - 1/n^2, {n, 2, Infinity}]
Out[1]= 0.5

In[2]:= NProduct[(n^2)/(n^2 - 1), {n, 2, Infinity}]
Out[2]= 2.0

In[3]:= NProduct[1 + 1/n^2, {n, 1, Infinity}]
Out[3]= 3.67608
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
