# FindSpanningTree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindSpanningTree[g] gives a spanning tree (forest) of g as a graph.`**

## Examples

_No verified examples yet for this function._

## Algorithm

spanningtree.c - FindSpanningTree[g].

A BFS spanning forest of the underlying undirected graph: for each component, the tree edges chosen by BFS are collected in their original form (preserving DirectedEdge/UndirectedEdge and orientation). Returns Graph[verts, treeEdges]; for a connected graph the tree has VertexCount - 1 edges.

Memory (SPEC section 4): returns a freshly-allocated Graph; frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
