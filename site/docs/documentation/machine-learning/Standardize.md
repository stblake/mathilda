# Standardize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Standardize[data] shifts each column of data to zero mean and rescales it to unit sample standard deviation (divisor n-1, matching StandardDeviation). A flat list is treated as n observations of one variable. A constant column becomes exactly 0 rather than Indeterminate.`**

## Examples (3)

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

## Implementation notes

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
