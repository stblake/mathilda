---
source: src/stats/mediandeviation.c
---
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
