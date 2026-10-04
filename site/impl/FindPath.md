---
source: src/graph/gops_cycles.c
---
**Algorithm.** `builtin_find_path` takes `FindPath[g, s, t]`, `FindPath[g, s, t, kspec]` and `FindPath[g, s, t, kspec, n]`, with `kspec` as for `FindCycle` (length counts edges) and `n` a positive integer or `All`. The basic form is an iterative DFS from `s` (`dfs_path`) with visited marks and neighbours taken in `EdgeList` order, returning `{path}` as a vertex list, or `{}` when `t` is unreachable. `s == t` gives `{}`. Bounded or multi-path forms backtrack over simple paths with an on-path mark, record every path that reaches `t` within the length range, and then order them shortest first (ties by discovery order).

**Data structures.** A cached forward-arc CSR built from the edge-index views (undirected edges contribute both arcs); `int` arrays for the DFS stack of (vertex, cursor) pairs, the current path and the on-path marks. The result is a list of vertex lists.

**Complexity / limits.** `O(V + E)` for the single-path form, and a repeated query on one graph reuses the cached CSR. Path enumeration is exponential in the worst case and is bounded by a `TimeConstrained`-aware step budget, after which the call stays unevaluated. A vertex not in the graph, or a weighted graph with a `kspec`, is also left unevaluated.
