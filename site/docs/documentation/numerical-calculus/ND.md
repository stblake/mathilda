# ND

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ND[expr, x, x0]`**

gives a numerical approximation to the derivative of expr with respect to x at the point x0.

**`ND[expr, {x, n}, x0]`**

gives a numerical approximation to the n-th derivative.

**`ND[{e1, e2, ...}, x, x0]`**

threads element-wise over the first argument.

<details>
<summary>Notes</summary>

Default Method -\> EulerSum uses Richardson extrapolation of forward, direction-Scale finite differences (works for non-analytic expr; needs integer n \>= 0). Method -\> NIntegrate uses Cauchy's integral formula via NResidue (needs expr analytic near x0; allows fractional/complex order). ND cannot recognize small numbers that should be zero -- Chop if needed. Options: Method (EulerSum | NIntegrate), Scale (step size / contour radius / complex direction, default 1), Terms (EulerSum starting extrapolation depth, default 7; grown adaptively to meet AccuracyGoal), WorkingPrecision, AccuracyGoal (default MachinePrecision), PrecisionGoal, MaxRecursion.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= ND[Exp[x], x, 1]
Out[1]= 2.71828

In[2]:= ND[Cos[x]^3, {x, 2}, 0]
Out[2]= -3.0

In[3]:= ND[Sin[x], x, Pi I]
Out[3]= 11.592 + 2.15751e-13*I

In[4]:= ND[{Exp[x], Sin[x]}, x, 1]
Out[4]= {2.71828, 0.540302}
```

Non-analytic: use EulerSum

```mathematica
In[5]:= ND[Re[Cos[I y]], y, 1]
Out[5]= 1.1752
```

### Options (3)

```mathematica
In[6]:= ND[Abs[x], {x, 1}, 0, Scale -> 1 + I]
Out[6]= 0.707107 - 0.707107*I

In[7]:= ND[Sin[100 x], x, 0, Scale -> 1/100]
Out[7]= 100.0

In[8]:= ND[Exp[x^2], {x, 4}, 0, Method -> NIntegrate]
Out[8]= 12.0 - 9.99201e-16*I
```

### Applications (4)

```mathematica
In[9]:= ND[Sin[x], x, 1]
Out[9]= 0.540302

In[10]:= ND[Gamma[x], x, 1]
Out[10]= -0.577216

In[11]:= ND[BesselJ[0, x], x, 2]
Out[11]= -0.576725

In[12]:= ND[Tan[x], {x, 2}, 1]
Out[12]= 11.4484
```

## Algorithm

nderiv.c — ND[expr, x, x0] / ND[expr, {x, n}, x0, opts]

Numerical approximation to the (n-th) derivative of `expr` w.r.t. `x` at x = x0. Two methods, selected by Method (default EulerSum):

```text
  Method -> EulerSum   Richardson (Romberg/Neville) extrapolation of the
                       n-th *forward* finite difference taken along the
                       complex direction Scale:

                         D(h) = (1/(s h)^n) sum_{k=0}^n (-1)^{n-k} C(n,k)
                                             f(x0 + k s h),   s = Scale,
                         h_i = 2^-i  (i = 0 .. Terms-1),
                         T(i,0) = D(h_i),
                         T(i,j) = T(i,j-1)
                                  + (T(i,j-1) - T(i-1,j-1)) / (2^j - 1),
                         result = T(Terms-1, Terms-1).

                       The forward (one-sided) stencil along `s` is what
                       gives directional/one-sided derivatives — e.g. the
                       left/right derivatives of Abs and the complex
                       direction Scale -> 1 + I. The error expansion of a
                       forward difference runs in *all* powers of h, hence
                       the (2^j - 1) Richardson denominator (not 4^j - 1).
                       Works for non-analytic f (only samples along `s`).
                       Requires integer order n >= 1.

  Method -> NIntegrate Cauchy integral formula via the existing NResidue:
                         f^(n)(x0) = n! Res_{z=x0} f(z)/(z-x0)^(n+1)
                                   = Gamma(n+1) *
                                     NResidue[expr/(x-x0)^(n+1), {x, x0},
                                              Radius -> Scale].
                       Scale is the contour radius (default 1; NResidue's
                       own tiny 1/100 default would cause heavy 1/r^n
                       cancellation, so we always pass Radius -> Scale).
                       Gamma(n+1) (rather than n!) gives fractional /
                       complex order. Requires expr analytic near x0.
```

Options: Method, Scale (default 1; step for EulerSum / contour radius for NIntegrate; may be complex), Terms (default 7; EulerSum tableau depth), WorkingPrecision (MachinePrecision | digits -> MPFR), PrecisionGoal, MaxRecursion (NIntegrate only). Threads element-wise over a List in arg 1.

Memory: receives `res` owned by the evaluator; returns a fresh Expr* on success or NULL (unevaluated). Never frees `res`. Any temporary OwnValue created for the EulerSum sampler is removed on every return path.

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [NResidue](../../numerical-calculus/NResidue/), [Chop](../../elementary-functions/Chop/), [Abs](../../arithmetic/Abs/), [NIntegrate](../../numerical-calculus/NIntegrate/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [PrecisionGoal](../../other-advanced/PrecisionGoal/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_nderiv.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nderiv.c)

## Notes & additional examples

### Notes

`ND[expr, x, x0]` numerically differentiates `expr` at `x = x0`. The first case
recovers `Cos[1] = 0.540302`. The Gamma example gives `Gamma'[1] = -EulerGamma`,
since `PolyGamma[0, 1] = -EulerGamma`. The Bessel example uses the identity
`BesselJ[0, x]' = -BesselJ[1, x]`, so the value is `-BesselJ[1, 2]`. The
`{x, 2}` form takes the second derivative. The default `Method -> EulerSum`
applies Richardson extrapolation to finite differences; `Method -> NIntegrate`
uses Cauchy's integral formula and allows fractional or complex orders. `ND`
cannot recognise small numbers that should be zero — `Chop` if needed.
