# WeightedAdjacencyMatrix

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`WeightedAdjacencyMatrix[g] gives the adjacency matrix of g with each entry the corresponding edge's weight (0 where there is no edge). Equal to AdjacencyMatrix[g] when g has no EdgeWeight.`**

## Examples (6)

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

### Applications (2)

Weights fill the adjacency cells

```mathematica
In[5]:= WeightedAdjacencyMatrix[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}, EdgeWeight -> {5, 7}]]
Out[5]= {{0, 5, 0}, {5, 0, 7}, {0, 7, 0}}
```

Unweighted: every present edge is 1

```mathematica
In[6]:= WeightedAdjacencyMatrix[CycleGraph[3]]
Out[6]= {{0, 1, 1}, {1, 0, 1}, {1, 1, 0}}
```

## Algorithm

wtadjmat.c - WeightedAdjacencyMatrix[g]: dense adjacency matrix filled with per-edge weights instead of a literal 1.

Same algorithm as AdjacencyMatrix (adjmat.c): a DirectedEdge[a,b] sets M[a][b] = weight(a,b); an UndirectedEdge sets both M[a][b] and M[b][a]. Any entry with no edge is 0. For a graph with no EdgeWeight, every weight defaults to 1 (graph_resolve_edge_weights), so WeightedAdjacencyMatrix[g] == AdjacencyMatrix[g] exactly for an unweighted g.

Memory (SPEC section 4): returns a freshly-allocated matrix; frees res.

## Implementation notes

**Algorithm.** `builtin_weighted_adjacency_matrix` builds the dense `n x n`
adjacency matrix filled with per-edge **weights** instead of a literal `1`. It is
the same construction as `AdjacencyMatrix`: a `DirectedEdge[a, b]` sets
`M[a][b] = weight(a, b)`, an `UndirectedEdge` sets both `M[a][b]` and `M[b][a]`,
and any cell with no edge is `0`. The weights come from
`graph_resolve_edge_weights`, the shared resolver that returns `g`'s `EdgeWeight`
list or `{1, ..., 1}` when there is none — so for an unweighted graph
`WeightedAdjacencyMatrix[g]` equals `AdjacencyMatrix[g]` exactly.

**Data structures.** A flat `n*n` grid of `Expr*` cells, a `GraphVIdx` mapping
each vertex to its row/column index, and the resolved weight list. The grid is
assembled into a `List` of row `List`s; empty cells become `expr_new_integer(0)`.

**Complexity / limits.** `O(n^2)` to materialize the dense matrix plus `O(E)` to
place the weights. Weights are copied verbatim, so symbolic weights appear as-is
in the matrix. Returns `NULL` (unevaluated) for a non-graph argument.

- `Protected`. Equal to `AdjacencyMatrix[g]` exactly when `g` has no
  `EdgeWeight` (every weight defaults to `1`). Undirected edges put their weight
  in both symmetric positions.
- Mathematica returns a `SparseArray`; Mathilda returns a dense matrix.
- Unevaluated on a non-graph.

**Attributes:** `Protected`.

## References

**See also:** [EdgeWeight](../../graphs/EdgeWeight/), [SparseArray](../../data-structures/SparseArray/)

- Source: [`src/graph/wtadjmat.c`](https://github.com/stblake/mathilda/blob/main/src/graph/wtadjmat.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)

## Notes & additional examples

### Notes

`WeightedAdjacencyMatrix[g]` is the dense adjacency matrix with each present
edge's `EdgeWeight` in place of a `1`; absent edges are `0`. An undirected edge
fills both symmetric cells, a directed edge only `M[tail][head]`.

For a graph with no weights every edge defaults to weight `1`, so the result
coincides with `AdjacencyMatrix[g]`. The row and column order is `VertexList`
order.
