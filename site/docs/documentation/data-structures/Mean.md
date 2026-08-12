# Mean

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Mean[data] gives the mean estimate of the elements in data.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Mean[<|"a" -> 2, "b" -> 4, "c" -> 6|>]
Out[1]= 4
```

### Applications (4)

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

## Performance

Measured on arm64 Darwin at commit `2dea9cc05`.

| case | n | time |
|---|---:|---:|
| list of machine reals | 1,000 | 5 us |
| list of machine reals | 10,000 | 7 us |
| list of machine reals | 100,000 | 21 us |

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Quartiles over 2x10^6 | 17.3 s | 17.2 s | 17.4 s |
| MovingAverage window 100 | 17.2 s | 2.02 s | 4.33 s |
| Median over 2x10^6 | 10.2 s | 7.79 s | 13.3 s |
| Skewness over 2x10^6 | 0.618 s | 0.582 s | 3.52 s |
| Kurtosis over 2x10^6 | 0.572 s | 0.505 s | 3.21 s |
| StandardDeviation over 2x10^6 | 0.325 s | 0.284 s | 0.927 s |

## Implementation notes

**Algorithm.** `builtin_mean` first probes its argument with `MatrixQ`; if true it computes column-wise means via `apply_columnwise` (which is `Map[Mean, Transpose[matrix]]`). Otherwise it requires a `List` (`ListQ`). For a vector of length `n` it dispatches on element kinds: if any element is `EXPR_REAL`, it sums to a `double` and returns `expr_new_real(sum/n)`; if all elements are exact integers/rationals it accumulates the sum in `int64_t` numerator/denominator pairs (reducing by `gcd` each step) and returns `make_rational(sum_n, sum_d * n)`. Anything symbolic falls back to `(1/n) * (Plus @@ data)` built as `Times`/`Apply` nodes and re-evaluated.

**Limits.** The exact-rational accumulator uses fixed `int64_t` arithmetic, so it can overflow for large/many rationals (no GMP promotion in this path). Empty list returns `NULL`.

**Attributes:** `Protected`.

## See also

[Total](../../arithmetic/Total/), [Min](../../data-structures/Min/), [Max](../../data-structures/Max/)

## References

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_ndarray_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_reduce.c)

## Notes & additional examples

### Notes

`Mean[data]` is the arithmetic mean — the sum of the elements divided by their
count. It works symbolically as well as numerically: `Mean[{a, b, c}]` returns
the exact closed form `(a + b + c)/3`. Numeric data stays in exact rational
arithmetic, so `Mean[{1, 2, 3, 4}]` is `5/2` (not `2.5`) and the mean of the
first ten squares is `77/2`, with no round-off. Combined with generators like
`Table` and `Range`, `Mean` gives exact averages of structured data sets.
