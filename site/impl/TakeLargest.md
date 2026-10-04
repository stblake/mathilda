---
source: src/sort.c
---
**Algorithm.** `builtin_take_largest[list, n]` returns the `n` largest elements in
descending order. For a plain list it calls the shared `take_extreme(coll, f =
NULL, n, largest = true)`: it forms `(key, payload)` pairs with the key being the
element itself, `qsort`s them ascending by key (`sortby_pair_cmp`, canonical
order), and copies `k = min(n, len)` payloads walking down from the top. Over an
association the subjects are the values and the result is the association of the
matching entries.

When the first argument is a packed `NDArray`, it takes the `ndstruct_take_extreme`
fast path instead — an `O(n log k)` bounded-heap selection over the raw buffer, with
no boxed `Expr` per element — and degrades to delist-and-re-evaluate only if that
declines.

**Data structures.** `SortByPair` array of `n` `(key, payload)` `Expr` pairs for
the list path (freed after the result is built); a size-`k` heap over machine
numbers for the NDArray path. The input argument array is borrowed.

**Complexity / limits.** List path `O(n log n)` (a full sort). NDArray path
`O(n log k)` through the heap. `n` must be an explicit integer; if it exceeds the
length, all elements are returned.
