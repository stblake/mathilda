---
source: src/stats/meandeviation.c
---
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
