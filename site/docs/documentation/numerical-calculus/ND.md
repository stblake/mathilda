# ND

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ND[expr, x, x0]
    gives a numerical approximation to the derivative of expr with respect to x at the point x0.
ND[expr, {x, n}, x0]
    gives a numerical approximation to the n-th derivative.
ND[{e1, e2, ...}, x, x0]
    threads element-wise over the first argument.

Default Method -> EulerSum uses Richardson extrapolation of forward, direction-Scale finite differences (works for non-analytic expr; needs integer n >= 0). Method -> NIntegrate uses Cauchy's integral formula via NResidue (needs expr analytic near x0; allows fractional/complex order). ND cannot recognize small numbers that should be zero -- Chop if needed.

Options: Method (EulerSum | NIntegrate), Scale (step size / contour radius / complex direction, default 1), Terms (EulerSum extrapolation terms, default 7), WorkingPrecision, PrecisionGoal, MaxRecursion.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ND[Exp[x], x, 1]
Out[1]= 2.71828

In[2]:= ND[Cos[x]^3, {x, 2}, 0]
Out[2]= -3.0

In[3]:= ND[Sin[x], x, Pi I]
Out[3]= 11.592 + 1.32527e-10*I

In[4]:= ND[{Exp[x], Sin[x]}, x, 1]
Out[4]= {2.71828, 0.540302}

In[5]:= ND[Re[Cos[I y]], y, 1]          (* non-analytic: use EulerSum *)
Out[5]= 1.1752

In[6]:= ND[Abs[x], {x, 1}, 0, Scale -> 1 + I]
Out[6]= 0.707107 - 0.707107*I

In[7]:= ND[Sin[100 x], x, 0, Scale -> 1/100]
Out[7]= 100.0

In[8]:= ND[Exp[x^2], {x, 4}, 0, Method -> NIntegrate]
Out[8]= 12.0 - 3.3723e-15*I
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
In[1]:= ND[Sin[x], x, 1]
Out[1]= 0.540302
```

```mathematica
In[1]:= ND[Gamma[x], x, 1]
Out[1]= -0.577216
```

```mathematica
In[1]:= ND[BesselJ[0, x], x, 2]
Out[1]= -0.576725
```

```mathematica
In[1]:= ND[Tan[x], {x, 2}, 1]
Out[1]= 11.4484
```

### Notes

`ND[expr, x, x0]` numerically differentiates `expr` at `x = x0`. The first case
recovers `Cos[1] = 0.540302`. The Gamma example gives `Gamma'[1] = -EulerGamma`,
since `PolyGamma[0, 1] = -EulerGamma`. The Bessel example uses the identity
`BesselJ[0, x]' = -BesselJ[1, x]`, so the value is `-BesselJ[1, 2]`. The
`{x, 2}` form takes the second derivative. The default `Method -> EulerSum`
applies Richardson extrapolation to finite differences; `Method -> NIntegrate`
uses Cauchy's integral formula and allows fractional or complex orders. `ND`
cannot recognise small numbers that should be zero — `Chop` if needed.
