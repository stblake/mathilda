# NProduct

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NProduct[f, {i, imin, imax}]`**

gives a numerical approximation to the product of f for i from imin to imax.

**`NProduct[f, {i, imin, imax, di}] uses step di. imax may be Infinity. NProduct[f, {i, ...}, {j, ...}, ...] evaluates a multidimensional product (an inner bound may depend on an outer index). The index is localised (HoldAll). Evaluated as Exp[NSum[Log[f], ...]], so the NSum engine (Euler-Maclaurin for monotone factors, Wynn's epsilon otherwise) and its convergence test carry over. Machine or arbitrary precision via WorkingPrecision.`**

<details>
<summary>Notes</summary>

Options: Method (Automatic | EulerMaclaurin | WynnEpsilon), WorkingPrecision (default MachinePrecision), NProductFactors (leading factors taken explicitly, default 15), NProductExtraFactors, WynnDegree, VerifyConvergence (default True; a divergent product gives ComplexInfinity), AccuracyGoal, PrecisionGoal.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= NProduct[1 - 1/n^2, {n, 2, Infinity}]
Out[1]= 0.5

In[2]:= NProduct[(n^2)/(n^2 - 1), {n, 2, Infinity}]
Out[2]= 2.0

In[3]:= NProduct[1 + 1/n^2, {n, 1, Infinity}]
Out[3]= 3.67608
```

### Applications (2)

```mathematica
In[1]:= NProduct[1 - 1/n^2, {n, 2, Infinity}]
Out[1]= 0.5
```

```mathematica
In[1]:= NProduct[Cos[1/n], {n, 1, Infinity}]
Out[1]= 0.388536
```

## Algorithm

```text
nprod.c — NProduct[f, {i, imin, imax (, di)}, opts]   (see nprod.h)
```

Strategy -------- Per Keiper 1992 ("The N functions of Mathematica", G.15), a numerical product is evaluated as the exponential of a numerical sum of logarithms:

```text
    Prod_{i=imin}^{imax} f(i)  =  Exp[ NSum[ Log[f(i)], {i, imin, imax} ] ].
```

This is exact for finite ranges (Exp inverts Log on any branch and the principal-log phases that wind by multiples of 2*pi are unwrapped by Exp) and correct for convergent infinite products (factors -> 1 => Log f -> 0 is smooth on the principal branch, and Euler-Maclaurin uses only branch-independent

```text
derivatives f'/f).  We therefore delegate every hard part — method selection,
```

Euler-Maclaurin, Wynn epsilon, Cohen-Villegas-Zagier, MPFR working precision, large-finite tail differences, and divergence detection — to the existing, tested NSum engine, and only:

```text
  - parse NProduct's own option names and map them onto NSum's;
  - add guard digits on the arbitrary-precision path (Exp amplifies the
    absolute error of the exponent into relative error of the product);
  - handle multidimensional products by recursion (inner NProduct as body);
  - special-case divergence (NSum -> ComplexInfinity).
```

Memory: receives `res` owned by the evaluator; returns a fresh Expr* on

```text
success or NULL (unevaluated).  Never frees `res`.
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## See also

[AccuracyGoal](../../other-advanced/AccuracyGoal/), [PrecisionGoal](../../other-advanced/PrecisionGoal/), [Exp](../../elementary-functions/Exp/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_nprod.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nprod.c)

## Notes & additional examples

### Notes

`NProduct[f, {i, imin, imax}]` numerically evaluates a product, with `imax`
allowed to be `Infinity`. The first case is a telescoping product:
`Product[1 - 1/n^2, {n, 2, Infinity}] = 1/2` exactly. The second converges to
about `0.388536`. `NProduct` is evaluated internally as `Exp[NSum[Log[f], ...]]`,
so the Euler–Maclaurin / Wynn's-epsilon machinery and convergence test of `NSum`
carry over. With `VerifyConvergence -> True` (default) a divergent product gives
`ComplexInfinity`. Use `WorkingPrecision` for arbitrary precision.
