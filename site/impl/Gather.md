---
references:
  - "T. H. Cormen, C. E. Leiserson, R. L. Rivest and C. Stein, *Introduction to Algorithms*, 3rd ed. (MIT Press, 2009), §11 (hash tables, open addressing)."
source: src/list/gather.c
---
**Algorithm.** `builtin_gather` (in `src/list/gather.c`) partitions a list into
sublists of structurally identical elements: two elements share a sublist iff
`expr_eq` holds, sublists appear in order of first occurrence, and input order is
kept within each sublist. Unlike `Split`, grouping is not restricted to adjacent
runs (`Gather[{1, 7, 3, 7, 2, 3, 9}]` → `{{1}, {7, 7}, {3, 3}, {2}, {9}}`). The
grouping is not reimplemented here — it calls `assoc_gather_core` (`assoc.c`)
with a `NULL` key function, which selects the identity key.

**Data structures.** The shared hash-indexed engine's `KeyIndex`
(open-addressing over `Expr*`) plus per-group doubling buffers; the identity path
takes a copy of each element as its own group key rather than evaluating
`Identity[x]`.

**Complexity / limits.** `O(n)` amortised. Passing `NULL` rather than
`f = Identity` makes `Gather[l] === GatherBy[l, Identity]` structural rather than
coincidental and avoids `n` redundant `Identity[x]` evaluations. Requires a
single argument; returns `NULL` otherwise.
