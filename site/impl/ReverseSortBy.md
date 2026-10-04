---
source: src/sort.c
---
**Algorithm.** `builtin_reverse_sort_by` is `Reverse` of `SortBy`: it re-heads its
arguments onto `SortBy` (`sort_call_as`), calls `builtin_sort_by`, and flips the
top level with the shared `reverse_top_level` helper. The three-argument form
`ReverseSortBy[coll, f, p]` is treated as `SortBy[coll, f, p[#2, #1] &]`
(`sort_by_with_p` with the ordering reversed), which differs from
`Reverse[SortBy[coll, f, p]]` on ties — again the Mathematica-15 convention.

As with `ReverseSort`, the `NDArray` arm of `reverse_top_level` is a correctness
fix: a packed buffer is reversed by a row-sized `memcpy` swap rather than silently
passing through as `SortBy`'s ascending result. Over an association the sort is by
`f` of each value.

**Data structures.** A synthesised `SortBy[...]` call (copies of the arguments)
and an in-place reversal of the sorted result's argument array (or a row-block
`memcpy` swap for a packed buffer).

**Complexity / limits.** `O(n)` key evaluations plus the `O(n log n)` underlying
`SortBy` and an `O(n)` reversal.
