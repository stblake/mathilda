# MedianDeviation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MedianDeviation[data]`**

gives the median absolute deviation from the median of the elements in data.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= MedianDeviation[{1, 2, 3, 4}]
Out[1]= 1

In[2]:= MedianDeviation[{1, 2, 3, 10}]
Out[2]= 1
```

## Algorithm

mediandeviation.c -- MedianDeviation[].

MedianDeviation[data] = Median[Abs[data - Median[data]]] (the median absolute deviation), composed through the evaluator so exact input stays exact. Matrix input recurses columnwise; NDArray input is materialised via pack_unpack (correctness-first, no kernel yet). See stats.h and stats_common.h for the subsystem layout.

## Implementation notes

- `Protected`.
- Exact input gives exact output: `MedianDeviation[{1, 2, 3, 4}]` is `1`.
- For `MatrixQ` data the deviation is computed per column.
- Requires numeric data (`MedianDeviation::rectn` otherwise).

**Attributes:** `Protected`.

## References

**See also:** [MatrixQ](../../expression-information/MatrixQ/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_quantile_family.c`](https://github.com/stblake/mathilda/blob/main/tests/test_quantile_family.c)
