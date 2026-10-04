# SeriesCoefficient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SeriesCoefficient[f, {x, x0, k}]`**

gives the coefficient of (x - x0)^k in the power-series expansion of f about x = x0. Works for a concrete integer index k and a finite expansion

<details>
<summary>Notes</summary>

point, for any f that Series can expand. Protected; the expansion variable must evaluate to a symbol.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= SeriesCoefficient[BesselJ[0, x], {x, 0, 4}]
Out[1]= 1/64

In[2]:= SeriesCoefficient[Exp[x], {x, 0, 5}]
Out[2]= 1/120
```

### Applications (5)

The coefficient of x^10 in e^x is 1/10!

```mathematica
In[3]:= SeriesCoefficient[Exp[x], {x, 0, 10}]
Out[3]= 1/3628800
```

Odd coefficients of Tan give the tangent numbers

```mathematica
In[4]:= SeriesCoefficient[Tan[x], {x, 0, 7}]
Out[4]= 17/315
```

The Fibonacci generating function: F(10) = 89

```mathematica
In[5]:= SeriesCoefficient[1/(1 - x - x^2), {x, 0, 10}]
Out[5]= 89
```

Even coefficients of Cos are reciprocals of factorials

```mathematica
In[6]:= SeriesCoefficient[Cos[x], {x, 0, 8}]
Out[6]= 1/40320
```

A symbolic index returns the closed-form general term

```mathematica
In[7]:= SeriesCoefficient[ProductLog[x], {x, 0, n}]
Out[7]= Piecewise[{{(-n)^(-1 + n)/Factorial[n], n >= 1}}, 0]
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

**Algorithm.** `builtin_seriescoefficient` implements `SeriesCoefficient[f, {x, x0, k}]` —
the coefficient of `(x − x0)^k` — by reusing `Series`. It first type-checks the expansion
variable with `series_spec_vars_ok` (the same gate `Series` uses: the variable must evaluate
to a symbol), placed *before* the symbolic-index special cases so that, with e.g. `x = 5` in
scope, a spec `{5, 0, n}` cannot spuriously match. Two closed-form general terms are then
recognised directly for a symbolic index `n` at `x0 = 0` and emitted as a `Piecewise`:
`ProductLog[x]` gives `(−n)^(n−1)/n!` for `n ≥ 1`, and `FresnelC`/`FresnelS[x]` give their
`4m+1` / `4m+3` power-class coefficients. Inexact arguments route through
`internal_rationalize_then_numericalize`. The general path parses the spec with
`parse_series_spec`, declines the leading `{x, x0}` two-argument form, expands with
`Series[f, {x, x0, max(k, 0)}]` (`eval_and_free`), and reads the coefficient out of the
resulting `SeriesData`: with minimal power `nmin` and denominator `den`, power `k` sits at
list index `j = k·den − nmin`, so an in-range `j` copies that coefficient and an out-of-range
`j` returns `0`. If `Series` collapses to something free of `x`, the result is that value for
`k = 0` and `0` otherwise.

**Data structures.** `Expr` trees throughout. The intermediate is the `SeriesData` six-tuple
`{var, x0, {c_0, …}, nmin, nmax, den}`; the index arithmetic (`k·den − nmin`) is done in
`long long`. The two symbolic-index cases build a `Piecewise[{{value, cond}}, 0]` `Expr`.
`SeriesCoefficient` is registered `Protected` and, unlike Mathematica folklore, does **not**
hold its arguments — the variable is type-checked instead (see the `series_init` comment),
so a named spec and a `Sequence @@` splice both work.

**Complexity / limits.** Cost is that of one `Series` expansion to order `k` (a full Taylor
pass through the series engine), followed by `O(1)` extraction; the result is exact (rational
or symbolic). Composite results that are a prefactor times a `SeriesData` (e.g. asymptotic
expansions at `Infinity`, `1/x` expansions), and non-integer indices, are left unevaluated
(`NULL`). The symbolic-index general term is produced only for the three heads above; any
other symbolic `n` declines. Being a structural head returning symbolic coefficients, it
carries no NDArray/packed kernel and no `Compile[]` lowering.

- `Protected` only; the arguments are evaluated and the expansion variable must be
  a symbol (see `Series` above — the same `ivar` decline applies, under the head
  `SeriesCoefficient`).
- Computed by expanding with `Series` and extracting the `k`-th coefficient from
  the resulting `SeriesData`; general for any head `Series` can expand, with a
  concrete integer index `k` and a finite expansion point.
- Composite results (a prefactor times a `SeriesData`, e.g. asymptotic
  expansions at Infinity) and non-integer indices are left unevaluated; the
  symbolic-index general term (a Piecewise) is not produced.

**Attributes:** `Protected`.

## References

**See also:** [Series](../../power-series/Series/), [SeriesData](../../power-series/SeriesData/)

- Source: [`src/calculus/series.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/series.c)
- Specification: [`docs/spec/builtins/power-series.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/power-series.md)
- Tests: [`tests/test_besselj.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besselj.c)
- Tests: [`tests/test_fresnelc.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fresnelc.c)
- Tests: [`tests/test_fresnels.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fresnels.c)
- Tests: [`tests/test_productlog.c`](https://github.com/stblake/mathilda/blob/main/tests/test_productlog.c)

## Notes & additional examples

### Notes

`SeriesCoefficient[f, {x, x0, k}]` returns the coefficient of `(x - x0)^k` in the
power-series expansion of `f` about `x = x0`, for any `f` that `Series` can expand.
At a concrete integer index it expands `f` to order `k` and extracts the single
coefficient, so the result is exact (rational or symbolic). For a handful of heads
(`ProductLog`, `FresnelC`, `FresnelS`) a **symbolic** index returns the closed-form
general term as a `Piecewise`. The arguments are evaluated normally
(`SeriesCoefficient` is `Protected`, not `HoldAll`); a non-integer or otherwise
unusable index is left unevaluated.
