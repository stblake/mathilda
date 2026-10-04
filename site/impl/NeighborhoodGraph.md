---
source: src/graph/gops_edit.c
---
**Algorithm.** `builtin_neighborhood_graph` returns the subgraph of `g` induced by
a set of centres together with every vertex within graph distance `k` of them
(default `k = 1`; `k = Infinity` and any non-negative integer are accepted, edge
**direction ignored**). The centre argument is a single vertex or a list of them;
non-vertices are silently dropped. From each centre it runs a BFS over the
all-direction incidence (`gops_inc_build` with `GOPS_INC_ALL`), stamped per centre
so the balls of several centres do not interfere, stopping when a vertex's distance
reaches `k`. Result vertices are listed centres-first, then each centre's newly
reached ball in `VertexList` order — a `qsort` of the ball, or a linear rescan over
all vertices when the ball exceeds roughly `1/8` of them. `build_induced` then keeps
every edge of `g` whose both endpoints were selected.

**Data structures.** A `GopsView` exposes the memo's endpoint arrays (`eu`/`ev`/
`edir`); `gops_inc_build` turns them into a CSR incidence (`GopsInc`, neighbours of
each vertex in one `start`/`nbr` pair). The BFS uses per-centre `mark` (a stamp),
`dist`, `queue`, and `ball` arrays plus a `pos[]` map from graph vertex to result
position. `build_induced` emits each kept edge once — at its later endpoint, in
lower-triangular adjacency order — reusing the shared edge, weight, and vertex nodes
through `expr_copy`, and carries an `EdgeWeight` list aligned to the surviving edges.

**Complexity / limits.** `O(V + E)` per centre for the BFS and `O(V + E)` to induce
the subgraph; each edge is emitted once. `k = Infinity` yields the whole connected
hull of the centres. Arity other than two or three, or a `k` that is neither a
non-negative integer nor `Infinity`, returns `NULL`. `res` is borrowed and never
modified.
