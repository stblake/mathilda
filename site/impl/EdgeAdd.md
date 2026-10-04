---
source: src/graph/gops_edit.c
---
**Algorithm.** `builtin_edge_add` accepts `EdgeAdd[g, e]` or `EdgeAdd[g, {e1, e2, ...}]`, where each item must parse as an edge (`gops_parse_edge`: `DirectedEdge`, `UndirectedEdge`, or the `->`/`<->` sugar); anything else leaves the call unevaluated. It copies the vertex, edge and endpoint arrays with room for the new edges, then walks the new edges in order. An endpoint already in the graph is found through the memoized vertex index; an unseen endpoint is appended to the vertex list once, in first-appearance order, via a small scratch hash. `u -> v` sugar takes the graph's kind: undirected when the graph has no directed edges, directed otherwise, while an explicit `DirectedEdge` stays directed. A weighted graph keeps its `EdgeWeight` list aligned, and each new edge gets weight 1.

**Data structures.** The graph is an ordinary `Graph[List[vertices], List[edges]]` expression tree with no dedicated `EXPR_*` tag. The edit works on the per-graph memo (`gops_view`: endpoint arrays `eu[k]`/`ev[k]`/`directed[k]`) and builds the result through `gops_graph_new`, sharing existing vertex and edge nodes and rebuilding only sugared edges.

**Complexity / limits.** `O(V + E)` for the copy plus one hash probe per new endpoint. Mathilda graphs are simple, so an edit that would create a self-loop or a parallel edge, such as adding an edge that already exists, is rejected when the result is validated and the call stays unevaluated instead of returning a multigraph.
