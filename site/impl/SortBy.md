---
source: src/sort.c
---
**Algorithm.** `builtin_sort_by` orders a collection by `f` applied to each element
(by `f[value]` for an association, the keys following). `f` is evaluated once per
element and the resulting key is cached in a `SortByPair`, then `qsort`
(`sortby_pair_cmp`) compares by that key via `expr_compare`, breaking ties by the
element itself and finally by original position — a stable canonical order. Two
variants extend this: a *list* key `{f1, f2, …}` builds the tuple `{f1[e], …}` and
sorts lexicographically; the three-argument `SortBy[coll, f, p]` first puts elements
in canonical order of their subjects, then merge-sorts by the ordering function `p`
on the `f`-values. The one-argument `SortBy[f]` returns the operator form
`Function[SortBy[#, f]]`.

**Data structures.** A `SortByPair{key, payload, pos, multi}` array; the C library
`qsort` for the two-argument form, a merge sort (`p_sort_perm`) for the `p` form.

**Complexity / limits.** O(n log n) comparisons with exactly one `f` evaluation per
element (the key is computed up front, not re-evaluated inside the comparator). The
`{f1, f2, …}` form sorts by `f1`, then `f2`, … as tie-breakers.
