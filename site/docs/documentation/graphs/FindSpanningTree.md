# FindSpanningTree

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindSpanningTree[g] gives a spanning tree (forest) of g as a graph.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= EdgeList[FindSpanningTree[CycleGraph[4]]]
Out[1]= {1 <-> 2, 4 <-> 1, 2 <-> 3}

In[2]:= EdgeCount[FindSpanningTree[CompleteGraph[6]]]
Out[2]= 5

In[3]:= EdgeList[FindSpanningTree[Graph[{1,2,3,4},{1->2,2->3,3->1}]]]
Out[3]= {1 -> 2, 3 -> 1}

In[4]:= EdgeList[FindSpanningTree[Graph[{1,2,3,4},{1<->2,3<->4}]]]
Out[4]= {1 <-> 2, 3 <-> 4}
```

## Algorithm

spanningtree.c - FindSpanningTree[g].

A BFS spanning forest of the underlying undirected graph: for each component, the tree edges chosen by BFS are collected in their original form (preserving DirectedEdge/UndirectedEdge and orientation). Returns Graph[verts, treeEdges]; for a connected graph the tree has VertexCount - 1 edges.

Memory (SPEC section 4): returns a freshly-allocated Graph; frees res.

## Implementation notes

- `Protected`. Has `VertexCount - 1` edges when `g` is connected; a
  disconnected graph gives a spanning forest. Tree edges keep their original
  direction. Unweighted (not a minimum-weight tree); see `FindShortestPath` for
  the shared machinery.
- Unevaluated on a non-graph.

**Attributes:** `Protected`.

## References

**See also:** [FindShortestPath](../../graphs/FindShortestPath/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
