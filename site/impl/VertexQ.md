---
source: src/graph/membership.c
---
**Algorithm.** `builtin_vertex_q` answers `VertexQ[g, v]` — `True` iff `v` is a
vertex of `g`. It calls `graph_vertex_position(g, v)` and returns the boolean
`position >= 0`. Membership is structural, matching the Wolfram Language's
`SameQ` semantics via `expr_eq`: the integer vertex `1` is *not* matched by the
real `1.0`. A `g` that is not a valid graph makes `graph_vertex_position` return
`-2`, so `VertexQ` gives `False` rather than leaving itself unevaluated. The head
requires exactly two arguments; any other arity returns `NULL`.

**Data structures.** The position lookup is an `O(1)` probe into the
validated-graph memo (`graph_util.c`), which already holds the `GraphVIdx`
`expr_hash` index of vertices that validating `g` built — no per-call scan of the
vertex list. The only allocation is the returned `True`/`False` symbol.

**Complexity / limits.** `O(1)` expected on a memo hit; the memo itself is built
in `O(V + E)` the first time a graph is validated. Because equality is
structural, callers that want numeric-insensitive membership must normalise
their vertices first.
