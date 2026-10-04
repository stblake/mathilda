# Standardize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Standardize[data] shifts each column of data to zero mean and rescales it to unit sample standard deviation (divisor n-1, matching StandardDeviation). A flat list is treated as n observations of one variable. A constant column becomes exactly 0 rather than Indeterminate.`**

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Standardize[{1., 2., 3., 4.}]
Out[1]= {-1.1619, -0.387298, 0.387298, 1.1619}

In[2]:= Standardize[{{1., 10.}, {2., 20.}, {3., 30.}}]
Out[2]= {{-1.0, -1.0}, {0.0, 0.0}, {1.0, 1.0}}

In[3]:= Standardize[{{1., 5.}, {2., 5.}, {3., 5.}}]
Out[3]= {{-1.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}}
```

### Applications (3)

A flat list is n observations of one variable

```mathematica
In[4]:= Standardize[{1., 2., 3., 4., 5.}]
Out[4]= {-1.26491, -0.632456, 0.0, 0.632456, 1.26491}
```

Each column is standardised independently

```mathematica
In[5]:= Standardize[{{1., 10.}, {2., 20.}, {3., 30.}}]
Out[5]= {{-1.0, -1.0}, {0.0, 0.0}, {1.0, 1.0}}
```

A constant column becomes exactly 0, not Indeterminate

```mathematica
In[6]:= Standardize[{{1., 5.}, {2., 5.}, {3., 5.}}]
Out[6]= {{-1.0, 0.0}, {0.0, 0.0}, {1.0, 0.0}}
```

## Implementation notes

**Algorithm.** `builtin_standardize` reads its argument into a row-major `n × dim`
buffer with `ml_read_data` — a flat list loads as `n` observations of *one* variable
(`was_vector`), a matrix as `n` rows of `dim` columns — and calls `ml_standardize(x, n,
dim, rescale = true)`. That kernel makes two passes over the columns: `ml_column_mean`
accumulates each column mean, then `ml_column_sd` accumulates the **sample** standard
deviation with the `n − 1` divisor (the same divisor `StandardDeviation` uses, so the
result agrees with `(x − Mean[x])/StandardDeviation[x]` written by hand). Each entry is
then shifted by its column mean and divided by its column standard deviation. A column
of zero variance is the one special case: rather than divide by zero and propagate a
`NaN` through every downstream row reduction, the kernel leaves such a column at exactly
`0.0` (`v = (sd[j] > 0.0) ? v / sd[j] : 0.0`).

**Data structures.** One contiguous `double` buffer from `ml_read_data`
(`na_load_matrix` / `na_load_vector`), plus a `dim`-length mean and SD scratch. The
result is rebuilt as a plain `List` through `ml_list_of_reals` / `ml_list_matrix`, so
the evaluator's own packing gate — not the ML code — decides whether it is held as a
buffer.

**Complexity / limits.** `O(n · dim)` time, two passes, `O(dim)` extra space. The mean
buffer is `calloc`-ed only to settle a phantom-uninitialised-read warning on the
degenerate `dim == 0` path; it carries no runtime cost.

- **Columns are variables, rows are observations.** A flat list is treated as `n`
  observations of *one* variable, not one observation of `n`.
- The divisor is `n - 1` (the sample standard deviation), matching
  `StandardDeviation` — so `Standardize[x]` agrees with
  `(x - Mean[x])/StandardDeviation[x]` written out by hand. A mismatch here would be
  invisible on the mean but not on the scale.
- **A constant column becomes exactly `0`, not `Indeterminate`.** Zero variance
  carries no information, so "no deviation from the mean" is the honest value;
  dividing by the zero standard deviation would propagate `Indeterminate` through
  every reduction over the row.

**Attributes:** `Protected`.

## References

**See also:** [StandardDeviation](../../data-structures/StandardDeviation/)

- Source: [`src/ml/pca.c`](https://github.com/stblake/mathilda/blob/main/src/ml/pca.c)
- Specification: [`docs/spec/builtins/machine-learning.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/machine-learning.md)
- Tests: [`tests/test_ml_pca.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ml_pca.c)

## Notes & additional examples

### Notes

Columns are variables and rows are observations: a flat list is `n` observations of a
*single* variable, not one observation of `n`.

The divisor is `n − 1`, the sample standard deviation, so `Standardize[x]` agrees to the
last digit with `(x − Mean[x])/StandardDeviation[x]` written out by hand. A constant
column carries no information, so it is shifted to exactly `0` rather than divided by its
zero standard deviation — which would propagate `Indeterminate` through every reduction
over the row.
