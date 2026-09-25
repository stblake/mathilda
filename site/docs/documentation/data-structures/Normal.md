# Normal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Normal[expr]`**

converts expr to a normal expression. If expr is a SeriesData object, the O-term is dropped and the truncated polynomial (or Laurent/Puiseux sum) is returned. Other expressions pass through unchanged.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Normal[{<|"a" -> 1|>, f[<|"b" -> 2|>]}]
Out[1]= {{"a" -> 1}, f[{"b" -> 2}]}

In[2]:= Normal[<|"a" -> <|"b" -> 1|>|>]
Out[2]= {"a" -> <|"b" -> 1|>}
```

### Applications (4)

```mathematica
In[3]:= Normal[Series[Exp[x], {x, 0, 5}]]
Out[3]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + 1/120 x^5

In[4]:= Normal[Series[Sin[x]/x, {x, 0, 6}]]
Out[4]= 1 - 1/6 x^2 + 1/120 x^4 - 1/5040 x^6

In[5]:= Normal[Series[Tan[x], {x, 0, 7}]]
Out[5]= x + 1/3 x^3 + 2/15 x^5 + 17/315 x^7

In[6]:= Normal[Series[Log[1 + x], {x, 0, 5}]]
Out[6]= x - 1/2 x^2 + 1/3 x^3 - 1/4 x^4 + 1/5 x^5
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

## Implementation notes

**Algorithm.** `builtin_normal` converts a `SeriesData[x, x0, {a0,...,a_{k-1}},
nmin, nmax, den]` into an ordinary polynomial by dropping the O-term. It builds
the base `(x - x0)` (`series_build_xmx0`), then for each non-zero coefficient `a_i`
forms the term `a_i (x - x0)^((nmin+i)/den)` — using an integer exponent when
`den == 1`, otherwise `Rational[num, den]`, and emitting the coefficient bare when
the exponent is 0 — and sums the terms (`Plus`, or the single term / literal `0`
for degenerate cases), evaluating the result. Any argument that is not a 6-element
`SeriesData` is passed through unchanged (`expr_copy`).

**Data structures.** A direct read of the `SeriesData` arg slots (coefficient
`List`, `nmin`, `den`); no `SeriesObj` is reconstructed.

**Attributes:** `Protected`.

## References

**See also:** [Hold](../../expression-information/Hold/), [HoldComplete](../../expression-information/HoldComplete/)

- Source: [`src/calculus/series.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/series.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_bessely.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bessely.c)
- Tests: [`tests/test_blas.c`](https://github.com/stblake/mathilda/blob/main/tests/test_blas.c)

## Notes & additional examples

### Notes

`Normal[expr]` converts `expr` to a normal expression. Applied to a `SeriesData`
object it drops the `O`-term and returns the truncated polynomial (or
Laurent/Puiseux sum) as an ordinary `Plus` expression, ready to be added,
differentiated, or substituted into. Expressions that are already normal pass
through unchanged.
