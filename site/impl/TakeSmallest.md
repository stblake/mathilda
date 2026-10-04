---
source: src/sort.c
---
**Algorithm.** `builtin_take_smallest[list, n]` returns the `n` smallest elements
in ascending order — the mirror of `TakeLargest`. For a plain list it calls the
shared `take_extreme(coll, f = NULL, n, largest = false)`: `(key, payload)` pairs
keyed by the element, `qsort`ed ascending (canonical order), then `k = min(n,
len)` payloads copied from the bottom. Over an association the subjects are the
values and the result is the association of the matching entries.

A packed `NDArray` first argument takes the `ndstruct_take_extreme` fast path
(`largest = false`), an `O(n log k)` bounded-heap selection over the raw buffer,
degrading to delist-and-re-evaluate only if that declines.

**Data structures.** `SortByPair` array of `(key, payload)` `Expr` pairs for the
list path; a size-`k` heap over machine numbers for the NDArray path. The input
argument array is borrowed.

**Complexity / limits.** List path `O(n log n)`; NDArray path `O(n log k)`. `n`
must be an explicit integer; a larger `n` returns every element.
