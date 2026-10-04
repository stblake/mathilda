---
source: src/graph/generators.c
---
**Algorithm.** `builtin_path_graph` builds an undirected path. `PathGraph[n]` with an integer
argument gives the vertices `1..n` with edges `i <-> i+1` for `i = 1..n-1`. `PathGraph[{v1, ...,
vk}]` with a list argument uses the given expressions as vertices (copied) and joins consecutive
ones, `v_i <-> v_{i+1}`. The function assembles a raw `Graph[List[verts], List[edges]]`
expression, which the evaluator then re-runs through `builtin_graph` to canonicalise and
validate.

**Data structures.** The vertex and edge `List`s are built directly; all edges are
`UndirectedEdge` nodes. No weights.

**Complexity / limits.** `O(n)`. The integer form needs a non-negative integer (negative or
non-integer returns unevaluated); `PathGraph[0]` is the empty path.
