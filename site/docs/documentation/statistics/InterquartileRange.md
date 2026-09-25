# InterquartileRange

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`InterquartileRange[data]`**

gives the difference between the upper and lower quartiles of the elements in data.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= InterquartileRange[{1, 2, 3, 4, 5, 6, 7, 8}]
Out[1]= 4

In[2]:= InterquartileRange[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}}]
Out[2]= {6, 6, 6}
```

## Algorithm

interquartilerange.c -- InterquartileRange[].

InterquartileRange[data] = q3 - q1 of Quartiles[data] (Wolfram's IQR uses the Quartiles parameterization {{1/2,0},{0,1}}, so composing over the Quartiles builtin is the definition, not a shortcut). Matrix input recurses per-column BEFORE the vector path: a k-column matrix's Quartiles result is a k-list of triples, and for k == 3 a bare "is it a 3-list" guard would confuse it with a vector's {q1,q2,q3} and compute Quartiles(col3)-Quartiles(col1) -- the silent wrong answer the STATS-1 plan review flagged. The vector-path guard also requires all three quartiles to be SCALARS for the same reason. NDArray input is materialised via pack_unpack (correctness-first, no kernel). See stats.h and stats_common.h for the subsystem layout.

## Implementation notes

- `Protected`.
- A robust scale estimator: insensitive to outliers beyond the quartiles.
- For `MatrixQ` data the IQR is computed per column.
- Requires numeric data (`InterquartileRange::rectn` otherwise).

**Attributes:** `Protected`.

## References

**See also:** [Quartiles](../../statistics/Quartiles/), [MatrixQ](../../expression-information/MatrixQ/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_quantile_family.c`](https://github.com/stblake/mathilda/blob/main/tests/test_quantile_family.c)
