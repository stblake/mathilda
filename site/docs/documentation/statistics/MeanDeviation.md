# MeanDeviation

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MeanDeviation[data]`**

gives the mean absolute deviation from the mean of the elements in data.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= MeanDeviation[{1, 2, 3, 4}]
Out[1]= 1

In[2]:= MeanDeviation[{1/2, 3/2}]
Out[2]= 1/2
```

### Applications (4)

The exact mean absolute deviation of an integer sample

```mathematica
In[3]:= MeanDeviation[{1, 2, 3, 4}]
Out[3]= 1
```

Exact rationals stay exact

```mathematica
In[4]:= MeanDeviation[{1/2, 3/2}]
Out[4]= 1/2
```

A machine-real sample gives a machine-real deviation

```mathematica
In[5]:= MeanDeviation[{2., 4., 4., 4., 5., 5., 7., 9.}]
Out[5]= 1.5
```

On a matrix the deviation is computed per column

```mathematica
In[6]:= MeanDeviation[{{1, 2}, {3, 4}, {5, 9}}]
Out[6]= {4/3, 8/3}
```

## Algorithm

meandeviation.c -- MeanDeviation[].

MeanDeviation[data] = Mean[Abs[data - Mean[data]]], composed through the evaluator so exact input stays exact (the src/stats exactness discipline -- see mean.c's overflow note). Matrix input recurses columnwise; NDArray input is materialised via pack_unpack (correctness-first, no kernel yet). See stats.h and stats_common.h for the subsystem layout.

## Implementation notes

**Algorithm.** `MeanDeviation[data]` is the mean absolute deviation from the mean,
`Mean[Abs[data − Mean[data]]]`. `builtin_meandeviation` builds exactly that one
tree and hands it to the evaluator, so the src/stats exactness discipline holds:
exact input gives exact output. The steps:

1. materialise a visible `NDArray` / packed argument via `pack_unpack`
   (correctness-first, no kernel yet);
2. require `data` to be a non-empty `List`;
3. matrix input (first element a `List`) recurses columnwise via
   `stats_apply_columnwise("MeanDeviation", data)`;
4. validate every element — it must be real-numeric **and** finite
   (`has_nonfinite` rejects `Infinity`, `ComplexInfinity`, `Indeterminate`, which
   are `NumericQ` and free of the imaginary unit and so would pass the element gate
   but then leave the composed `Mean` tree half-evaluated); `MeanDeviation::rectn`
   otherwise;
5. build and evaluate `Mean[Abs[data − Mean[data]]]` — the inner `Mean`, the
   `Listable` `Plus`/`Abs`, and the outer `Mean` all reduce because the data was
   just verified all-numeric.

**Data structures.** A single composed `Expr` tree; no auxiliary buffers. The
matrix recursion reuses `stats_apply_columnwise`
(`Map[MeanDeviation, Transpose[matrix]]`).

**Complexity / limits.** `O(n)` to build, then the evaluator's two `Mean`
reductions over the vector. No packed buffer kernel and no `Compile[]` lowering —
a visible `NDArray` is unpacked to the exact `List` path first.

- `Protected`.
- Exact input gives exact output: `MeanDeviation[{1, 2, 3, 4}]` is `1`.
- For `MatrixQ` data the deviation is computed per column.
- Requires numeric data (`MeanDeviation::rectn` otherwise).

**Attributes:** `Protected`.

## References

**See also:** [MatrixQ](../../expression-information/MatrixQ/)

- Source: [`src/stats/meandeviation.c`](https://github.com/stblake/mathilda/blob/main/src/stats/meandeviation.c)
- Specification: [`docs/spec/builtins/statistics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/statistics.md)
- Tests: [`tests/test_quantile_family.c`](https://github.com/stblake/mathilda/blob/main/tests/test_quantile_family.c)

## Notes & additional examples

### Notes

`MeanDeviation[data] = Mean[Abs[data − Mean[data]]]`, the mean absolute deviation
about the mean. It is composed through the evaluator, so exact input gives exact
output. Non-finite entries (`Infinity`, `ComplexInfinity`, `Indeterminate`) are
rejected even though they are `NumericQ` and free of the imaginary unit — they
would otherwise leave the composed `Mean` half-evaluated — and non-numeric data
triggers `MeanDeviation::rectn`.

Matrix input is reduced per column; a visible `NDArray` is materialised to the
exact `List` path.
