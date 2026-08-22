# Median

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Median[data]`**

gives the median estimate of the elements in data.

**`Median[dist]`**

gives the median of the distribution dist.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Median[<|"a" -> 1, "b" -> 3, "c" -> 5|>]
Out[1]= 3

In[2]:= Variance[<|"a" -> 2, "b" -> 4, "c" -> 6|>]
Out[2]= 4

In[3]:= StandardDeviation[<|"a" -> 2, "b" -> 4, "c" -> 6|>]
Out[3]= 2
```

### Applications (3)

```mathematica
In[4]:= Median[{5, 1, 3, 2, 4}]
Out[4]= 3

In[5]:= Median[{1, 2, 3, 4}]
Out[5]= 5/2

In[6]:= Median[Table[k^2, {k, 1, 10}]]
Out[6]= 61/2
```

## Performance

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

**Algorithm.** `builtin_median` requires a `List`. If the first element is itself a `List` it treats the input as a matrix/tensor and reduces column-wise through `apply_columnwise` (`Map[Median, Transpose[...]]`). For a 1-D vector it first verifies every element is a real numeric via the helper `is_real_numeric` (which checks `NumericQ` and `FreeQ[#, I]`); non-real data prints `Median::rectn` and leaves the call unevaluated. It then evaluates `Sort[data]`: for odd `n` it returns the middle element (`sorted[n/2]`); for even `n` it returns `(sorted[n/2-1] + sorted[n/2]) / 2`, built as `Plus` then `Divide` and re-evaluated so the result stays exact (rational) when the inputs are exact. `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## References

**See also:** [Variance](../../data-structures/Variance/), [StandardDeviation](../../data-structures/StandardDeviation/), [Mean](../../data-structures/Mean/)

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_bignum_rational_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bignum_rational_numeric.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_ndarray_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_reduce.c)

## Notes & additional examples

### Notes

`Median[data]` returns the middle value of the sorted data. For an odd number of
elements it is the single central element (the first example sorts to
`{1, 2, 3, 4, 5}`, giving `3`); for an even number it is the exact average of the
two central elements, so `Median[{1, 2, 3, 4}]` is `5/2`. The result is kept in
exact arithmetic — the median of the first ten squares is `61/2`, the average of
the 5th and 6th sorted values `25` and `36`. `Median` expects numeric data; an
even-length list of unresolved symbols cannot be averaged and is left
unevaluated with a `Median::rectn` message.
