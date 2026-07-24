# NLimit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NLimit[expr, z -> z0]
    numerically finds the limiting value of expr as z approaches z0.

A geometric sequence of sample points approaching z0 is constructed (z0 may be finite, complex, or an infinite point such as Infinity or I Infinity) and the limit is recovered by sequence acceleration. Method -> Automatic (default) runs both Richardson/Romberg and Wynn's epsilon and keeps the most self-consistent estimate, so branch-point / fractional-power approaches (which Richardson alone mishandles) are resolved accurately (Levin's u-transform also participates when the samples are settling). Method -> EulerSum forces Richardson/Romberg; Method -> SequenceLimit forces Wynn's epsilon; Method -> "Levin" forces Levin's transformation ("LevinU" | "LevinT" | "LevinV" select the u/t/v variant). expr must be numerical when z is numerical. Small spurious residuals are not recognised as zero -- Chop if needed.

Options: Method (Automatic | EulerSum | SequenceLimit | "Levin"), WorkingPrecision (default MachinePrecision), Direction (Automatic == -1, or a complex approach vector), Scale (initial step / distance, default 1), Terms (default 7), WynnDegree (SequenceLimit iterations, default 1).
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

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= NLimit[Sin[x]/x, x -> 0]
Out[1]= 1.0
```

```mathematica
In[1]:= NLimit[(1 + 1/n)^n, n -> Infinity]
Out[1]= 2.71828
```

```mathematica
In[1]:= NLimit[n (2^(1/n) - 1), n -> Infinity]
Out[1]= 0.693147
```

```mathematica
In[1]:= NLimit[Zeta[x] - 1/(x - 1), x -> 1]
Out[1]= 0.577216
```

```mathematica
In[1]:= NLimit[Sin[x]/x, x -> 0, Method -> "Levin"]
Out[1]= 1.0
```

### Notes

`NLimit[expr, z -> z0]` builds a geometric sequence of sample points approaching
`z0` and recovers the limit by sequence acceleration. The first three cases give
`1`, the constant `E = 2.71828...`, and `Log[2] = 0.693147...`. The fourth is the
classic Laurent-expansion limit of the Riemann zeta function at its pole: the
constant term is the Euler–Mascheroni constant `EulerGamma = 0.577216...`. `z0`
may be finite, complex, or an infinite point such as `Infinity` or `I Infinity`.
`Method -> Automatic` (default) keeps the most self-consistent of Richardson
extrapolation (`EulerSum`), Wynn's epsilon (`SequenceLimit`) and Levin's
u-transform; `Method -> "Levin"` (`"LevinU"`/`"LevinT"`/`"LevinV"`) forces
Levin's transformation. `Chop` small spurious residuals.
