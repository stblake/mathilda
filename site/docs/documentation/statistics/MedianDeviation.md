# MedianDeviation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MedianDeviation[data]`**

gives the median absolute deviation from the median of the elements in data.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= MedianDeviation[{1, 2, 3, 4}]
Out[1]= 1

In[2]:= MedianDeviation[{1, 2, 3, 10}]
Out[2]= 1
```

### Applications (3)

The exact median absolute deviation (MAD) of an integer sample

```mathematica
In[3]:= MedianDeviation[{1, 2, 3, 4}]
Out[3]= 1
```

A machine-real sample gives a machine-real MAD

```mathematica
In[4]:= MedianDeviation[{1., 1., 2., 2., 4., 6., 9.}]
Out[4]= 1.0
```

On a matrix the MAD is computed per column

```mathematica
In[5]:= MedianDeviation[{{1, 2}, {3, 4}, {5, 9}}]
Out[5]= {2, 2}
```

## Algorithm

mediandeviation.c -- MedianDeviation[].

MedianDeviation[data] = Median[Abs[data - Median[data]]] (the median absolute deviation), composed through the evaluator so exact input stays exact. Matrix input recurses columnwise; NDArray input is materialised via pack_unpack (correctness-first, no kernel yet). See stats.h and stats_common.h for the subsystem layout.

## Implementation notes

**Algorithm.** `MedianDeviation[data]` is the median absolute deviation (MAD) from
the median, `Median[Abs[data − Median[data]]]` — a robust scale estimator.
`builtin_mediandeviation` mirrors `MeanDeviation` with the center and the outer
reduction both `Median` instead of `Mean`:

1. materialise a visible `NDArray` / packed argument via `pack_unpack`
   (correctness-first, no kernel yet);
2. require `data` to be a non-empty `List`;
3. matrix input (first element a `List`) recurses columnwise via
   `stats_apply_columnwise("MedianDeviation", data)`;
4. validate every element — real-numeric **and** finite; `has_nonfinite` rejects
   `Infinity` / `ComplexInfinity` / `Indeterminate` (which would pass the
   `NumericQ` + no-`I` element gate but leave the composed `Median` tree
   half-evaluated); `MedianDeviation::rectn` otherwise;
5. build and evaluate `Median[Abs[data − Median[data]]]`.

Composing through the evaluator keeps exact input exact — `MedianDeviation[{1, 2,
3, 4}]` is the exact `1`.

**Data structures.** A single composed `Expr` tree; no auxiliary buffers. The
matrix recursion reuses `stats_apply_columnwise`
(`Map[MedianDeviation, Transpose[matrix]]`).

**Complexity / limits.** Two `Median` reductions over the vector (each `O(n log
n)` through `Median`'s own sort). No packed buffer kernel and no `Compile[]`
lowering — a visible `NDArray` is unpacked to the exact `List` path first.

- `Protected`.
- Exact input gives exact output: `MedianDeviation[{1, 2, 3, 4}]` is `1`.
- For `MatrixQ` data the deviation is computed per column.
- Requires numeric data (`MedianDeviation::rectn` otherwise).

**Attributes:** `Protected`.

## References

**See also:** [MatrixQ](../../expression-information/MatrixQ/)

- Source: [`src/stats/mediandeviation.c`](https://github.com/stblake/mathilda/blob/main/src/stats/mediandeviation.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_quantile_family.c`](https://github.com/stblake/mathilda/blob/main/tests/test_quantile_family.c)

## Notes & additional examples

### Notes

`MedianDeviation[data] = Median[Abs[data − Median[data]]]`, the median absolute
deviation (MAD) about the median — a robust scale estimator whose breakdown point
is far higher than that of the standard deviation. Composed through the evaluator,
so exact input gives exact output.

As with `MeanDeviation`, non-finite entries (`Infinity`, `ComplexInfinity`,
`Indeterminate`) are rejected rather than left half-evaluated, non-numeric data
triggers `MedianDeviation::rectn`, matrix input is reduced per column, and a
visible `NDArray` is materialised to the exact `List` path.
