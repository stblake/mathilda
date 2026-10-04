---
source: src/graph/gmet_generators.c
---
**Algorithm.** `builtin_kary_tree` builds a `k`-ary tree on `n` vertices with the
vertices numbered in **breadth-first** order: `KaryTree[n]` is binary (`k = 2`),
`KaryTree[n, k]` is `k`-ary. The helper `kary_tree(n, k)` emits, for each vertex
`i` (0-based internally), the edges to its children `k*i + 1, ..., k*i + k` as
long as the child index is `< n`. This is the implicit-heap layout, so the parent
of vertex `c` is `(c - 1) / k`; the first vertices fill complete levels and the
last level is filled left to right. (The sibling head `CompleteKaryTree[n, k]`
builds the *complete* `k`-ary tree of `n` levels by summing the level sizes and
calling the same helper.)

**Data structures.** A growable `Pairs` buffer of integer endpoint pairs; the
finished graph is assembled by `pairs_graph(n, &p)`, which produces the canonical
`Graph` on vertices `1..n`, undirected, and seeds the memo. No adjacency or
hashing is needed — the tree structure is arithmetic on the indices.

**Complexity / limits.** `O(n)` — one edge per non-root vertex. `n` is bounded by
`GEN_MAX_VERTICES`; a negative `n` or `k`, or an over-cap size, leaves the call
unevaluated.
