# Kurtosis

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Kurtosis[data]`**

gives the coefficient of kurtosis (peak/tail vs flank concentration) of data, equivalent to CentralMoment\[data, 4\] / CentralMoment\[data, 2\]^2. For a matrix it is taken columnwise.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Kurtosis[{1, 2, 3, 4, 5}]
Out[1]= 17/10

In[2]:= Kurtosis[{1, 2, 4, 8}]
Out[2]= 25141/13225
```

## Algorithm

kurtosis.c -- Kurtosis[]. Split from stats.c; see stats.h and stats_common.h for the subsystem layout.

Kurtosis[data] -- the coefficient of kurtosis, a measure of peak/tail vs flank concentration. Equivalent to CentralMoment[data, 4] / CentralMoment[data, 2]^2 (Pearson kurtosis, not the excess form). For a matrix or array it is taken columnwise (the CentralMoment ratio threads). The shared body lives in stats_common.c (stats_standardized_moment), which also routes NDArray / packed inputs to the buffer kernel ndred_kurtosis.

## Implementation notes

- `Protected`.
- Equivalent to `CentralMoment[data, 4] / CentralMoment[data, 2]^2` (Pearson kurtosis, not the excess form).
- For a matrix, gives the columnwise kurtoses.
- Handles numerical and symbolic data; exact input gives exact output.
- Fast path on `NDArray`/packed real buffers (`ndred_kurtosis`) and lowerable inside `Compile[]`; an integer buffer degrades to the exact `List` result.

**Attributes:** `Protected`.

## See also

[NDArray](../../linear-algebra/NDArray/), [List](../../other-advanced/List/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_stats.c`](https://github.com/stblake/mathilda/blob/main/tests/test_stats.c)
