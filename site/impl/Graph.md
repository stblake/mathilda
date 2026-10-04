---
source: src/graph/construct.c
---
**Algorithm.** `builtin_graph` is the canonicalising constructor. A graph that is
already valid is its own fixed point, so `graph_is_valid` is checked first and the
builtin returns `NULL` for it — skipping a full rebuild-and-compare on every
re-evaluation. Otherwise `try_build_canonical` tells the two input forms apart by
the second argument: `Graph[verts, edges, opt...]` when `args[1]` is a `List` (the
explicit edge list), `Graph[edges, opt...]` otherwise (`args[1]` is the first
option). Edge sugar is normalised on construction — `Rule`/`u -> v` becomes
`DirectedEdge[u, v]`, `TwoWayRule`/`u <-> v` becomes `UndirectedEdge[u, v]`, and an
already-normalised edge passes through. Per-edge options `EdgeWeight` and
`EdgeCapacity` are each a `Rule[key, List]` of length equal to the edge count,
slotted by canonical rank so every spelling of one graph canonicalises to the same
tree (`EdgeWeight` before `EdgeCapacity`, whatever order they were written).

**Data structures.** The object is a plain `Expr` tree
`Graph[List verts, List edges, opt...]` — no new `EXPR_*` tag. Vertices are copied
from an explicit list or derived from the normalised edges in first-appearance
order; in the same pass every endpoint is resolved to a vertex index through a
`GraphVIdx` open-addressing hash (the linear `expr_eq` rescan it replaces made
deriving the vertex list `O(E·V)`, seconds for a 40000-edge graph). The index and
the endpoint arrays (`eu`/`ev`/`edir`) are handed to `graph_memo_seed`, which
validates self-loops and parallel edges and memoises the result keyed on the node
pointer (kept alive and immutable), so the first `VertexQ`/`EdgeQ` on a freshly
built graph is an `O(1)` memo hit. One `DirectedEdge` and one `UndirectedEdge` head
node is cached per construction and shared by `expr_copy` into every rebuilt edge.

**Complexity / limits.** Construction is `O(V + E)` with one hash per endpoint.
Mathilda graphs are **simple**: a self-loop, a parallel edge, a 3-argument edge, an
endpoint absent from an explicit vertex list, or a wrong-length / repeated /
unknown option leaves `Graph[...]` unevaluated (`NULL`). The canonical tree is built
entirely from `expr_copy` of the argument's parts, so `res` is never cannibalised;
the already-canonical case returns `NULL` to reach a fixed point.
