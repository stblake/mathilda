---
source: src/sort.c
---
**Algorithm.** `builtin_take_smallest_by[list, f, n]` returns the `n` elements for
which `f` is smallest, in ascending order of `f` — the mirror of `TakeLargestBy`.
It calls the shared `take_extreme(coll, f, n, largest = false)`: a `(key,
payload)` pair per element keyed by the eagerly evaluated `f[subject]` (for an
association, `f` of each value), `qsort`ed ascending in canonical order, and the
bottom `k = min(n, len)` payloads copied from the smallest key up. The original
elements are returned, keeping the collection's head.

**Data structures.** `SortByPair` array of `n` `(key, payload)` `Expr` pairs; each
key is a freshly built and evaluated `f[subject]`. The pairs are freed after the
result is built; the input argument array is borrowed.

**Complexity / limits.** `O(n)` key evaluations plus an `O(n log n)` sort. `n`
must be an explicit integer; a larger `n` returns every element.
