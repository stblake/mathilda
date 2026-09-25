# MeanDeviation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MeanDeviation[data]`**

gives the mean absolute deviation from the mean of the elements in data.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= MeanDeviation[{1, 2, 3, 4}]
Out[1]= 1

In[2]:= MeanDeviation[{1/2, 3/2}]
Out[2]= 1/2
```

## Algorithm

meandeviation.c -- MeanDeviation[].

MeanDeviation[data] = Mean[Abs[data - Mean[data]]], composed through the evaluator so exact input stays exact (the src/stats exactness discipline -- see mean.c's overflow note). Matrix input recurses columnwise; NDArray input is materialised via pack_unpack (correctness-first, no kernel yet). See stats.h and stats_common.h for the subsystem layout.

## Implementation notes

- `Protected`.
- Exact input gives exact output: `MeanDeviation[{1, 2, 3, 4}]` is `1`.
- For `MatrixQ` data the deviation is computed per column.
- Requires numeric data (`MeanDeviation::rectn` otherwise).

**Attributes:** `Protected`.

## References

**See also:** [MatrixQ](../../expression-information/MatrixQ/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_quantile_family.c`](https://github.com/stblake/mathilda/blob/main/tests/test_quantile_family.c)
