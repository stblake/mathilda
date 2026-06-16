# SeriesCoefficient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SeriesCoefficient[f, {x, x0, k}]
    gives the coefficient of (x - x0)^k in the power-series expansion of f
    about x = x0. Works for a concrete integer index k and a finite expansion
point, for any f that Series can expand. HoldAll, Protected.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SeriesCoefficient[BesselJ[0, x], {x, 0, 4}]
Out[1]= 1/64

In[2]:= SeriesCoefficient[Exp[x], {x, 0, 5}]
Out[2]= 1/120
```

## Implementation notes

- `HoldAll`, `Protected`.
- Computed by expanding with `Series` and extracting the `k`-th coefficient from

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/power-series.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/power-series.md)

## Notes & additional examples

### Worked examples

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

### Notes

`SeriesCoefficient[f, {x, x0, k}]` returns the coefficient of `(x - x0)^k` in the
power-series expansion of `f` about `x = x0`, for any `f` that `Series` can expand
and a concrete integer index `k`. It is computed by expanding `f` to order `k` and
extracting the single coefficient, so the result is exact (rational or symbolic).
`SeriesCoefficient` is `HoldAll`, so the expansion variable is held unevaluated.
