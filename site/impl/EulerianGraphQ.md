---
source: src/graph/gops_preds.c
---
**Algorithm.** `builtin_eulerian_graph_q` (via `gops_eulerian`) is `True` when
`g` has a closed walk using every edge exactly once. The degree test is the
classical one: every vertex of even degree for an undirected graph, or
`in-degree = out-degree` for a directed graph, **and** all edges lying in one
connected component. One pass over the edges accumulates a per-vertex balance
(parity XOR for undirected, `+1`/`-1` for directed) and marks the touched
vertices; a union-find then checks that the non-isolated vertices form a single
component. For a directed graph, balanced plus weakly connected implies strongly
connected, so weak connectivity suffices. An edgeless graph with at least one
vertex is Eulerian; the null graph is not; a mixed graph is left unevaluated.

**Data structures.** A `GopsView` over the graph plus three `O(V)` scratch
arrays: the balance counters `bal[]`, a `touched[]` bitmap, and the union-find
`parent[]`. All work is on the view's integer endpoint arrays.

**Complexity / limits.** `O(V + E)`. A non-graph argument gives `False`; a mixed
graph returns `NULL` (unevaluated), as Mathematica leaves it. `FindEulerianCycle`
constructs an actual tour (Hierholzer) when this predicate holds.
