# Mean

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Mean[data] gives the mean estimate of the elements in data.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Mean[{1, 2, 3, 4}]
Out[1]= 5/2

In[2]:= Mean[{{a, u}, {b, v}, {c, w}}]
Out[2]= {1/3 (a + b + c), 1/3 (u + v + w)}
```

## Implementation notes

**Algorithm.** `builtin_mean` first probes its argument with `MatrixQ`; if true it computes column-wise means via `apply_columnwise` (which is `Map[Mean, Transpose[matrix]]`). Otherwise it requires a `List` (`ListQ`). For a vector of length `n` it dispatches on element kinds: if any element is `EXPR_REAL`, it sums to a `double` and returns `expr_new_real(sum/n)`; if all elements are exact integers/rationals it accumulates the sum in `int64_t` numerator/denominator pairs (reducing by `gcd` each step) and returns `make_rational(sum_n, sum_d * n)`. Anything symbolic falls back to `(1/n) * (Plus @@ data)` built as `Times`/`Apply` nodes and re-evaluated.

**Limits.** The exact-rational accumulator uses fixed `int64_t` arithmetic, so it can overflow for large/many rationals (no GMP promotion in this path). Empty list returns `NULL`.

- `Protected`.
- Supports numerical and symbolic data.
- For vectors, computes $(1/n) \sum x_i$.
- For matrices, computes means of elements in each column.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Mean[{1, 2, 3, 4}]
Out[1]= 5/2
```

```mathematica
In[1]:= Mean[{a, b, c}]
Out[1]= 1/3 (a + b + c)
```

```mathematica
In[1]:= Mean[{1/2, 1/3, 1/6}]
Out[1]= 1/3
```

```mathematica
In[1]:= Mean[Table[k^2, {k, 1, 10}]]
Out[1]= 77/2
```

### Notes

`Mean[data]` is the arithmetic mean — the sum of the elements divided by their
count. It works symbolically as well as numerically: `Mean[{a, b, c}]` returns
the exact closed form `(a + b + c)/3`. Numeric data stays in exact rational
arithmetic, so `Mean[{1, 2, 3, 4}]` is `5/2` (not `2.5`) and the mean of the
first ten squares is `77/2`, with no round-off. Combined with generators like
`Table` and `Range`, `Mean` gives exact averages of structured data sets.
