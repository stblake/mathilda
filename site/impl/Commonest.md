---
source: src/list.c
---
**Algorithm.** `builtin_commonest` (in `src/list.c`) tallies element multiplicities with a `HashTable` (one pass, O(N)), recording each distinct element's count and first-appearance index. The tallies are packed into `CommonestItem` records and sorted by descending count (ties broken by first occurrence) to find the maximum multiplicity. Without a count argument it returns every element sharing the maximum count; `Commonest[list, n]` / `Commonest[list, UpTo[n]]` returns up to the `n` most frequent.

**Data structures.** Parallel `Expr**` / `int64_t*` arrays for unique elements and their multiplicities, a `HashTable` mapping element → index, and a `CommonestItem` array (`element`, `count`, `first_index`) used for the stable sort.
