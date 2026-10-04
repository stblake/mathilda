---
source: src/graph/gops_transform.c
---
**Algorithm.** `builtin_graph_complement` gives the graph on the same vertices
with an edge wherever `g` has none. For an undirected `g` every non-adjacent pair
`i < j` becomes `i <-> j`; for a directed or mixed `g` every ordered pair `(i,
j)`, `i != j`, with no edge usable from `i` to `j` becomes the directed edge `i
-> j`. For each source vertex `i` it stamps `i`'s existing out-neighbours into a
scratch array and then emits an edge to every unstamped `j`, in row-major
`VertexList` order. Weights are dropped.

**Data structures.** A `GopsView` of the graph, an out-incidence CSR
(`gops_inc_build`), a per-row `stamp[]` array marking existing neighbours, and a
growable `EdgeBuf` sized up front to the exact complement edge count. Result
edges share a single interned edge head (`GopsHeads`) and the finished graph is
seeded into the memo (`gops_graph_new`).

**Complexity / limits.** Output-bound: `O(V^2)` since the complement of a sparse
graph is dense. The `tc_check_deadline` poll every 64 rows keeps it interruptible
under `TimeConstrained`.
