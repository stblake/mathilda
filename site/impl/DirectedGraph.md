---
source: src/graph/gops_transform.c
---
**Algorithm.** `builtin_directed_graph` converts a graph to a directed one. By default each
undirected edge `u <-> v` becomes the two directed edges `u -> v` and `v -> u` (weight
duplicated), while existing directed edges are kept in place; a fully directed graph is returned
unchanged. `DirectedGraph[g, "Acyclic"]` instead *orients* each undirected edge from the vertex
earlier in `VertexList` to the one later, yielding a DAG when the input is undirected (the edges
are then counting-sorted by tail then head); a mixed graph keeps its existing directed edges and
orders the rest.

**Data structures.** A `GopsView` with integer endpoints and a directed-edge count `ndir`; the
default form allocates `ne + (ne - ndir)` result edges, the "Acyclic" form reuses the edge count.
The result is a canonical `Graph` built by `gops_graph_new`.

**Complexity / limits.** `O(V + E)` (plus a counting sort in the undirected "Acyclic" case). A
second argument other than the string `"Acyclic"`, or a non-graph first argument, returns
unevaluated.
