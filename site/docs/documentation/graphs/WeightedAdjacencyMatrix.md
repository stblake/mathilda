# WeightedAdjacencyMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`WeightedAdjacencyMatrix[g] gives the adjacency matrix of g with each entry the corresponding edge's weight (0 where there is no edge). Equal to AdjacencyMatrix[g] when g has no EdgeWeight.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= WeightedAdjacencyMatrix[CycleGraph[4]] == AdjacencyMatrix[CycleGraph[4]]
Out[1]= True

In[2]:= WeightedAdjacencyMatrix[5]
Out[2]= WeightedAdjacencyMatrix[5]
```

### Options (2)

```mathematica
In[3]:= WeightedAdjacencyMatrix[Graph[{1,2,3},{1->2,2->3},EdgeWeight->{5,7}]]
Out[3]= {{0, 5, 0}, {0, 0, 7}, {0, 0, 0}}

In[4]:= WeightedAdjacencyMatrix[Graph[{1,2,3},{1<->2,2<->3},EdgeWeight->{5,7}]]
Out[4]= {{0, 5, 0}, {5, 0, 7}, {0, 7, 0}}
```

## Algorithm

wtadjmat.c - WeightedAdjacencyMatrix[g]: dense adjacency matrix filled with per-edge weights instead of a literal 1.

Same algorithm as AdjacencyMatrix (adjmat.c): a DirectedEdge[a,b] sets M[a][b] = weight(a,b); an UndirectedEdge sets both M[a][b] and M[b][a]. Any entry with no edge is 0. For a graph with no EdgeWeight, every weight defaults to 1 (graph_resolve_edge_weights), so WeightedAdjacencyMatrix[g] == AdjacencyMatrix[g] exactly for an unweighted g.

Memory (SPEC section 4): returns a freshly-allocated matrix; frees res.

## Implementation notes

- `Protected`. Equal to `AdjacencyMatrix[g]` exactly when `g` has no
  `EdgeWeight` (every weight defaults to `1`). Undirected edges put their weight
  in both symmetric positions.
- Mathematica returns a `SparseArray`; Mathilda returns a dense matrix.
- Unevaluated on a non-graph.

**Attributes:** `Protected`.

## References

**See also:** [EdgeWeight](../../graphs/EdgeWeight/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
