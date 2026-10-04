---
references:
  - "J. Bentley and B. Floyd, *Programming pearls: a sample of brilliance*, Communications of the ACM **30**(9) (1987) 754-757 (Floyd's sampling)."
source: src/graph/hyp_random.c
---
**Algorithm.** `builtin_random_hypergraph` has two families.
`RandomHypergraph[{n, m}, k]` (and the `..., c` repeat form) builds a `k`-uniform
hypergraph on `1..n`: each of the `m` hyperedges is an independent uniformly
random `k`-subset drawn by **Floyd's algorithm** — for `t = n-k .. n-1`, pick
`r` in `[0, t]` and take `r` unless it is already in the subset, in which case take
`t` — then sorted (insertion sort for `k <= 16`, else `qsort`). The
Function-Repository form `RandomHypergraph[{n, {e, a}}]` (or a list of `{e_i, a_i}`
pairs) draws vertices uniformly **with replacement**, so a vertex may repeat inside
a hyperedge; the vertex set is the vertices that actually occur, in first-appearance
order (the constructor derives them).

**Data structures.** A reusable `stamp` membership array (compared to the hyperedge
index `j`, so no clearing between draws) and a `buf` of chosen indices; the shared
`List` head and the `1..n` integer nodes are built once and shared by `expr_copy`.
All randomness is `random_uniform_01`, the user-visible stream, so both forms
reproduce under `SeedRandom`.

**Complexity / limits.** `O(k)` per `k`-uniform hyperedge (no rejection). Requires
`0 <= k <= n`, else unevaluated. Counts are capped at `INT32_MAX`. Where the FR
function returns the bare edge List, this returns a validated `Hypergraph`
(`EdgeList` recovers the List).
