# Median

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Median[data]
    gives the median estimate of the elements in data.
Median[dist]
    gives the median of the distribution dist.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Median[{1, 2, 3, 4, 5, 6, 7}]
Out[1]= 4

In[2]:= Median[{1, 2, 3, 4, 5, 6, 7, 8}]
Out[2]= 9/2

In[3]:= Median[{1, 2, 3, 4}]
Out[3]= 5/2

In[4]:= Median[{Pi, E, 2}]
Out[4]= E

In[5]:= Median[{1., 2., 3., 4.}]
Out[5]= 2.5

In[6]:= Median[{{1, 11, 3}, {4, 6, 7}}]
Out[6]= {5/2, 17/2, 5}

In[7]:= Median[{{{3, 7}, {2, 1}}, {{5, 19}, {12, 4}}}]
Out[7]= {{4, 13}, {7, 5/2}}

In[8]:= Median[{a, b, c}]
Out[8]= Median[{a, b, c}]
```

## Implementation notes

**Algorithm.** `builtin_median` requires a `List`. If the first element is itself a `List` it treats the input as a matrix/tensor and reduces column-wise through `apply_columnwise` (`Map[Median, Transpose[...]]`). For a 1-D vector it first verifies every element is a real numeric via the helper `is_real_numeric` (which checks `NumericQ` and `FreeQ[#, I]`); non-real data prints `Median::rectn` and leaves the call unevaluated. It then evaluates `Sort[data]`: for odd `n` it returns the middle element (`sorted[n/2]`); for even `n` it returns `(sorted[n/2-1] + sorted[n/2]) / 2`, built as `Plus` then `Divide` and re-evaluated so the result stays exact (rational) when the inputs are exact. `ATTR_PROTECTED`.

- `Protected`.
- Median is a robust location estimator, which means it not very sensitive to outliers.
- For `VectorQ` data $\{x_1, \dots, x_n\}$, the median can be thought of as the "middle value". Formally, when data is sorted as $\{x_{(1)}, \dots, x_{(n)}\}$, the median is given by the center element $x_{((n+1)/2)}$ if $n$ is odd and the mean of the two center elements $(x_{(n/2)} + x_{(n/2+1)})/2$ if $n$ is even.
- For `MatrixQ` data, the median is computed for each column vector. `Median` for a tensor gives columnwise medians at the first level.
- `Median` requires numeric values.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Median[{5, 1, 3, 2, 4}]
Out[1]= 3
```

```mathematica
In[1]:= Median[{1, 2, 3, 4}]
Out[1]= 5/2
```

```mathematica
In[1]:= Median[Table[k^2, {k, 1, 10}]]
Out[1]= 61/2
```

### Notes

`Median[data]` returns the middle value of the sorted data. For an odd number of
elements it is the single central element (the first example sorts to
`{1, 2, 3, 4, 5}`, giving `3`); for an even number it is the exact average of the
two central elements, so `Median[{1, 2, 3, 4}]` is `5/2`. The result is kept in
exact arithmetic — the median of the first ten squares is `61/2`, the average of
the 5th and 6th sorted values `25` and `36`. `Median` expects numeric data; an
even-length list of unresolved symbols cannot be averaged and is left
unevaluated with a `Median::rectn` message.
