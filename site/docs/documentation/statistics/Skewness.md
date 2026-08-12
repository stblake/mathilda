# Skewness

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Skewness[data]`**

gives the coefficient of skewness (a measure of asymmetry) of data, equivalent to CentralMoment\[data, 3\] / CentralMoment\[data, 2\]^(3/2). For a matrix it is taken columnwise.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Skewness[{1, 2, 3, 10}]
Out[1]= 18/25 Sqrt[2]

In[2]:= Skewness[{1., 2., 3., 4., 5.}]
Out[2]= 0.0
```

## Algorithm

skewness.c -- Skewness[]. Split from stats.c; see stats.h and stats_common.h for the subsystem layout.

Skewness[data] -- the coefficient of skewness, a measure of asymmetry. Equivalent to CentralMoment[data, 3] / CentralMoment[data, 2]^(3/2). For a matrix or array it is taken columnwise (the CentralMoment ratio threads). The shared body lives in stats_common.c (stats_standardized_moment), which also routes NDArray / packed inputs to the buffer kernel ndred_skewness.

## Implementation notes

- `Protected`.
- Equivalent to `CentralMoment[data, 3] / CentralMoment[data, 2]^(3/2)`.
- For a matrix, gives the columnwise skewnesses.
- Handles numerical and symbolic data; exact input gives exact output (a radical in general).
- Fast path on `NDArray`/packed real buffers (`ndred_skewness`) and lowerable inside `Compile[]`; an integer buffer degrades to the exact `List` result.

**Attributes:** `Protected`.

## See also

[NDArray](../../linear-algebra/NDArray/), [List](../../other-advanced/List/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)
