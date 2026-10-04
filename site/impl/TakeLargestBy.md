---
source: src/sort.c
---
**Algorithm.** `builtin_take_largest_by[list, f, n]` returns the `n` elements for
which `f` is largest, in descending order of `f`. It calls the shared
`take_extreme(coll, f, n, largest = true)`: a `(key, payload)` pair is formed per
element with the key `f[subject]` evaluated eagerly (for an association, `f` of
each value), the pairs are `qsort`ed ascending by key in canonical order, and the
top `k = min(n, len)` payloads are copied walking down from the largest key. The
original elements — not the keys — are returned, keeping the collection's head.

**Data structures.** `SortByPair` array of `n` `(key, payload)` `Expr` pairs; each
key is a freshly built and evaluated `f[subject]`. All pairs are freed once the
result node is built; the input argument array is borrowed.

**Complexity / limits.** `O(n)` key evaluations plus an `O(n log n)` sort (the
`By` form does not take the bounded-heap NDArray fast path that keyless
`TakeLargest` uses). `n` must be an explicit integer; a larger `n` returns every
element.
