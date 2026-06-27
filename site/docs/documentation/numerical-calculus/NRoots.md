# NRoots

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NRoots[lhs == rhs, var]
    yields a disjunction of equations var==r1 || var==r2 || ... giving numerical approximations to the roots of the polynomial equation in var. Roots of multiplicity k appear as k identical equations; a single root yields a bare equation. Real and complex coefficients are handled at machine and arbitrary precision. Method -> Automatic uses the Aberth-Ehrlich simultaneous iteration; "CompanionMatrix" uses companion-matrix eigenvalues (real QR directly, complex via a real 2n embedding); "JenkinsTraub" uses the three-stage Jenkins-Traub algorithm.

Options: Method (Automatic | "Aberth" | "CompanionMatrix" | "JenkinsTraub"), PrecisionGoal (Automatic = machine; a digit count selects arbitrary precision), MaxIterations, StepMonitor.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NRoots[1 + 2 x + 3 x^2 + 4 x^3 == 0, x]
Out[1]= x == -0.60583 || x == -0.0720852 - 0.638327*I || x == -0.0720852 + 0.638327*I

In[2]:= NRoots[x^2 - 2 == 0, x]
Out[2]= x == -1.41421 || x == 1.41421

In[3]:= NRoots[x^2 + 1 == 0, x]
Out[3]= x == 0.0 - 1.0*I || x == 0.0 + 1.0*I

In[4]:= NRoots[(x - 1)^3 == 0, x]
Out[4]= x == 1.0 || x == 1.0 || x == 1.0

In[5]:= NRoots[x^2 - (3 + 4 I) == 0, x]
Out[5]= x == -2.0 - 1.0*I || x == 2.0 + 1.0*I

In[6]:= NRoots[x^2 - 2 == 0, x, PrecisionGoal -> 30]
Out[6]= x == -1.414213562373095048801688724209 || x == 1.414213562373095048801688724209
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
