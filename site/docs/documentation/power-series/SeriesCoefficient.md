# SeriesCoefficient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SeriesCoefficient[f, {x, x0, k}]`**

gives the coefficient of (x - x0)^k in the power-series expansion of f about x = x0. Works for a concrete integer index k and a finite expansion

<details>
<summary>Notes</summary>

point, for any f that Series can expand. HoldAll, Protected.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= SeriesCoefficient[BesselJ[0, x], {x, 0, 4}]
Out[1]= 1/64

In[2]:= SeriesCoefficient[Exp[x], {x, 0, 5}]
Out[2]= 1/120
```

### Applications (4)

```mathematica
In[1]:= SeriesCoefficient[Exp[x], {x, 0, 10}]
Out[1]= 1/3628800
```

The coefficient of `x^7` in `Tan[x]` matches the corresponding tangent number:

```mathematica
In[1]:= SeriesCoefficient[Tan[x], {x, 0, 7}]
Out[1]= 17/315
```

The coefficient of `x^n` in `1/(1 - x - x^2)` is the n-th Fibonacci number; here
`F(10) = 89`:

```mathematica
In[1]:= SeriesCoefficient[1/(1 - x - x^2), {x, 0, 10}]
Out[1]= 89
```

```mathematica
In[1]:= SeriesCoefficient[Cos[x], {x, 0, 8}]
Out[1]= 1/40320
```

## Algorithm

============================================================================ series.c - Series and SeriesData ============================================================================

This module implements the power-series machinery for Mathilda.

SeriesData[x, x0, {a0, ..., a_{k-1}}, nmin, nmax, den] is the data head that represents a truncated power series. The i-th coefficient multiplies (x - x0)^((nmin + i)/den) and an O[x - x0]^(nmax/den) term captures the dropped higher-order terms.

```text
Series[f, {x, x0, n}]  expands f as a power series in (x - x0) up to
```

order n. Series also accepts the leading-term form Series[f, x -> x0] and the iterated multivariate form Series[f, {x, x0, nx}, {y, y0, ny}, ...]. The algorithm is a recursive "series algebra": primitive subexpressions become SeriesObj's, algebraic heads (Plus, Times, Power) combine them, and elementary heads (Exp, Log, Sin, Cos, Sinh, Cosh, Tan, Tanh) apply their known series kernels. Unknown heads fall back to naive Taylor via D[...]. Expansion about Infinity is handled by substituting x -> 1/u internally and presenting the result with Power[x, -1] as the series variable.

Normal[s] drops the O-term from a SeriesData and returns an ordinary sum.

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Limit x (Log[x+1]-Log[x]) at Infinity | 1.24 s | 2.09 s | -- |
| Limit Sin[x]/x at 0 | 0.963 s | 1.37 s | -- |
| Limit (Exp[x]-1-x)/x^2 at 0 | 0.218 s | 1.63 s | -- |
| Series Exp[Sin[x]] to order 20 | -- | -- | -- |
| Series 1/(1-x-x^2) to order 60 | -- | -- | -- |
| Series Log[1+Sin[x]] to order 24 | -- | -- | -- |

## Implementation notes

- `HoldAll`, `Protected`.
- Computed by expanding with `Series` and extracting the `k`-th coefficient from
  the resulting `SeriesData`; general for any head `Series` can expand, with a
  concrete integer index `k` and a finite expansion point.
- Composite results (a prefactor times a `SeriesData`, e.g. asymptotic
  expansions at Infinity) and non-integer indices are left unevaluated; the
  symbolic-index general term (a Piecewise) is not produced.

**Attributes:** `HoldAll`, `Protected`.

## See also

[HoldAll](../../expression-information/HoldAll/), [Series](../../power-series/Series/), [SeriesData](../../power-series/SeriesData/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/power-series.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/power-series.md)
- Tests: [`tests/test_besselj.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besselj.c)
- Tests: [`tests/test_fresnelc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fresnelc.c)
- Tests: [`tests/test_fresnels.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fresnels.c)
- Tests: [`tests/test_productlog.c`](https://github.com/stblake/mathilda/blob/main/tests/test_productlog.c)

## Notes & additional examples

### Notes

`SeriesCoefficient[f, {x, x0, k}]` returns the coefficient of `(x - x0)^k` in the
power-series expansion of `f` about `x = x0`, for any `f` that `Series` can expand
and a concrete integer index `k`. It is computed by expanding `f` to order `k` and
extracting the single coefficient, so the result is exact (rational or symbolic).
`SeriesCoefficient` is `HoldAll`, so the expansion variable is held unevaluated.
