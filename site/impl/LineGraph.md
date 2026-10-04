---
source: src/graph/gops_transform.c
---
**Algorithm.** `builtin_line_graph` builds the line graph: one vertex per edge of the input,
labelled by its `EdgeList` position `1..m` (rather than by the original edge expression). For an
undirected graph two line-graph vertices are joined when the corresponding edges share an
endpoint: scanning the incidence of each edge's endpoints, an undirected edge `j <-> i` is
emitted once, at the later edge. For a fully directed graph, edge `i -> edge j` whenever edge
`i` ends at the vertex where edge `j` starts. A mixed directed/undirected graph is left
unevaluated.

**Data structures.** A CSR incidence structure (`gops_inc_build`, modes `GOPS_INC_ALL` for
undirected and `GOPS_INC_OUT` for directed); the result is assembled as a canonical `Graph` on
integer vertices.

**Complexity / limits.** `O(E + sum of deg^2)` — the per-endpoint neighbour scans dominate, so a
high-degree vertex of degree `d` contributes about `d^2` line-graph edges. The result is a
canonical `Graph`; a non-graph or mixed argument returns unevaluated.
