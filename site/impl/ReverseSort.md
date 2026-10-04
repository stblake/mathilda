---
source: src/sort.c
---
**Algorithm.** `builtin_reverse_sort` is `Reverse` of `Sort`: it re-heads its
arguments onto `Sort` (`sort_call_as`, so the ascending routine runs as `Sort`
even if it re-evaluates itself), calls `builtin_sort`, and flips the top level with
`reverse_top_level`. The one exception is `ReverseSort[assoc]`, which sorts the
entries descending by value with **equal values kept in input order**
(`sort_collection_p`) — not `Reverse[Sort[assoc]]`, which would flip the ties,
matching Mathematica 15.

`reverse_top_level` handles both list representations. Its `NDArray` arm is a
correctness fix rather than an optimisation: an `NDArray` is `EXPR_NDARRAY`, not
`EXPR_FUNCTION`, so a plain early return left `ReverseSort[NDArray[...]]` showing
`Sort`'s ascending answer; it instead reverses whole rows in place via a
row-sized `memcpy` swap.

**Data structures.** A freshly built `Sort[...]` call (adopting copies of the
arguments), then an in-place pointer reversal of the sorted list's argument array,
or a row-block `memcpy` swap for a packed buffer.

**Complexity / limits.** `O(n log n)` from the underlying `Sort` plus an `O(n)`
reversal. Ordering is Mathilda's canonical order.
