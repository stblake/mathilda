# Normal

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Normal[expr]
    converts expr to a normal expression. If expr is a SeriesData object, the
    O-term is dropped and the truncated polynomial (or Laurent/Puiseux sum) is
    returned. Other expressions pass through unchanged.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Normal[Series[Exp[x], {x, 0, 5}]]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + 1/120 x^5

In[2]:= Normal[a + b]
Out[2]= a + b
```

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

- `Protected`.
- Returns the Plus of the coefficient-times-power terms (zero coefficients skipped). For non-`SeriesData` input, `Normal` is the identity.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/calculus/series.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/series.c)
- Specification: [`docs/spec/builtins/power-series.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/power-series.md)
