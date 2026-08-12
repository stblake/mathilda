# MovingMedian

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MovingMedian[list, r]`**

gives the moving median of list, computed using spans of r elements.

<details>
<summary>Notes</summary>

MovingMedian returns a list of length Length\[list\] - r + 1; for matrix input the medians are taken column-wise within each row-window. MovingMedian requires real-valued numeric data and stays unevaluated when r \< 1 or r \> Length\[list\].

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= MovingMedian[{1, 2, 5, 6, 1, 4, 3}, 3]
Out[1]= {2, 5, 5, 4, 3}

In[2]:= MovingMedian[{{1, 2}, {5, 3}, {1, 4}, {3, 2}, {5, 5}}, 2]
Out[2]= {{3, 5/2}, {3, 7/2}, {2, 3}, {4, 7/2}}

In[3]:= MovingMedian[N[{1, 5, 7, 3, 6, 2}], 3]
Out[3]= {5.0, 5.0, 6.0, 3.0}

In[4]:= MovingMedian[{1, 2, 3, 4}, 2]
Out[4]= {3/2, 5/2, 7/2}

In[5]:= MovingMedian[{2^100, 2^101, 2^102, 2^103}, 2]
Out[5]= {1901475900342344102245054808064, 3802951800684688204490109616128, 7605903601369376408980219232256}
```

### Applications (3)

```mathematica
In[1]:= MovingMedian[{1, 3, 2, 8, 5, 4, 9}, 3]
Out[1]= {2, 3, 5, 5, 5}
```

```mathematica
In[1]:= MovingMedian[{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, 5]
Out[1]= {3, 4, 5, 6, 7, 8}
```

```mathematica
In[1]:= MovingMedian[Table[Mod[n^2, 11], {n, 1, 12}], 4]
Out[1]= {9/2, 9/2, 4, 4, 4, 9/2, 9/2, 5/2, 1}
```

## Implementation notes

**Algorithm.** `builtin_moving_median` takes `(list, r)` with `r` a positive integer window (`EXPR_INTEGER`/`EXPR_BIGINT`); output length is `n - r + 1` and the call is unevaluated unless `1 <= r <= n`. It auto-detects vector vs. matrix mode by whether the first element is a `List`, and validates every leaf with `is_real_numeric` (matrices must be rectangular); invalid data prints `MovingMedian::arg1` and leaves the call unevaluated. It then slides the window, builds an `r`-element sublist, and delegates each window to `Median` (so a matrix window yields a column-wise median vector). Window results are assembled into the output `List`. `ATTR_PROTECTED`.

- `Protected`.
- Output length is `Length[list] - r + 1`.
- Operates on real-valued vectors and matrices. For matrix input, each window of `r` consecutive rows is reduced via `Median`, yielding a column-wise median vector per window.
- Exact rationals, bignums (arbitrary-precision integers), machine-precision reals, and `NumericQ`-real symbolic constants (`Pi`, `E`, ...) are all supported. Even-window medians yield exact rational midpoints when the data is exact.
- Stays unevaluated when `r < 1`, when `r > Length[list]`, when `r` is non-integer, or when the first argument is not a `List`.
- Non-numeric data triggers the `MovingMedian::arg1` message and the expression remains unevaluated.

**Attributes:** `Protected`.

## See also

[Median](../../data-structures/Median/), [NumericQ](../../expression-information/NumericQ/), [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/), [List](../../other-advanced/List/)

## References

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_ndarray_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_reduce.c)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)

## Notes & additional examples

### Notes

`MovingMedian[list, r]` returns a list of length `Length[list] - r + 1`: the
median of each window of `r` consecutive elements as the window slides across
`list`. Even-width windows average the two central order statistics, so results
may be exact rationals (e.g. `9/2`) rather than integers. The robust median
makes it less sensitive to outliers than `MovingAverage`.
