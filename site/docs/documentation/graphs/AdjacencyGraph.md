# AdjacencyGraph

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AdjacencyGraph[m] builds a graph on vertices 1..n from a 0/1 adjacency matrix m (undirected if m is symmetric, else directed).`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= AdjacencyGraph[{{0,1},{1,0}}]
Out[1]= Graph[<2 vertices, 1 edge>]

In[2]:= EdgeList[AdjacencyGraph[{{0,1,0},{0,0,1},{1,0,0}}]]
Out[2]= {1 -> 2, 2 -> 3, 3 -> 1}

In[3]:= EdgeList[AdjacencyGraph[{{0,1,1},{1,0,0},{1,0,0}}]]
Out[3]= {1 <-> 2, 1 <-> 3}

In[4]:= AdjacencyMatrix[AdjacencyGraph[AdjacencyMatrix[CycleGraph[4]]]] == AdjacencyMatrix[CycleGraph[4]]
Out[4]= True

In[5]:= AdjacencyGraph[{{1,2},{3}}]
Out[5]= AdjacencyGraph[{{1, 2}, {3}}]
```

## Algorithm

adjgraph.c - AdjacencyGraph[m]: build a graph from a 0/1 adjacency matrix.

The inverse of AdjacencyMatrix. Vertices are the integers 1..n. A symmetric matrix yields an undirected graph (one UndirectedEdge per i<j with m[i][j]=1); an asymmetric matrix yields a directed graph (a DirectedEdge for each off- diagonal m[i][j]=1). Diagonal entries (self-loops) are ignored. The result is returned as a Graph[...] expression and canonicalized/validated by the evaluator (builtin_graph).

Round-trips with AdjacencyMatrix when the source graph's vertices are 1..n.

Memory (SPEC section 4): returns a freshly-allocated Graph; frees res.

## Implementation notes

- `Protected`. The inverse of `AdjacencyMatrix`: undirected if `m` is symmetric,
  else directed. `AdjacencyGraph[AdjacencyMatrix[g]]` reproduces `g`'s edges
  (as a set — edges are regenerated in row-major order, so `EdgeList` order and
  undirected-edge orientation may differ from `g`'s; the adjacency matrices
  agree).
- A non-square (ragged) matrix, or entries other than 0/1, leave the call
  unevaluated. A `1` on the diagonal is currently dropped silently rather than
  rejected (the graph model has no self-loops).

**Attributes:** `Protected`.

## References

**See also:** [AdjacencyMatrix](../../graphs/AdjacencyMatrix/), [EdgeList](../../graphs/EdgeList/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
