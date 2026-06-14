# NLimit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NLimit[expr, z -> z0]
    numerically finds the limiting value of expr as z approaches z0.

A geometric sequence of sample points approaching z0 is constructed (z0 may be finite, complex, or an infinite point such as Infinity or I Infinity) and the limit is recovered by sequence acceleration. Method -> EulerSum (default) uses Richardson/Romberg extrapolation; Method -> SequenceLimit uses Wynn's epsilon algorithm. expr must be numerical when z is numerical. Small spurious residuals are not recognised as zero -- Chop if needed.

Options: Method (EulerSum | SequenceLimit), WorkingPrecision (default MachinePrecision), Direction (Automatic == -1, or a complex approach vector), Scale (initial step / distance, default 1), Terms (default 7), WynnDegree (SequenceLimit iterations, default 1).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NLimit[Sin[x]/x, x -> 0]
Out[1]= 1.0

In[2]:= NLimit[(1 + 1/n)^n, n -> Infinity]
Out[2]= 2.71828

In[3]:= NLimit[(1 + I/x)^x, x -> Infinity]
Out[3]= 0.540302 + 0.841471*I

In[4]:= NLimit[Tanh[Pi x]/(1 + x^2), x -> I] // Chop
Out[4]= 5.33707e-06 - 1.5708*I

In[5]:= NLimit[(10^x - 1)/x, x -> 0, Terms -> 10, Method -> SequenceLimit]
Out[5]= 2.30262

In[6]:= NLimit[z + Conjugate[z]/z, z -> 0, Direction -> -I] // Chop
Out[6]= -1.0

In[7]:= NLimit[Tan[z], z -> Infinity I, Method -> SequenceLimit] // Chop
Out[7]= 0.0 + 1.0*I

In[8]:= NLimit[(2^x - 1)/x, x -> 0, WorkingPrecision -> 30, Terms -> 14]
Out[8]= 0.6931471805599453094172321284473
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
