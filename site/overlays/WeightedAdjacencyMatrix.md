### Worked examples

```mathematica
In[1]:= WeightedAdjacencyMatrix[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}, EdgeWeight -> {5, 7}]]  (* weights fill the adjacency cells *)
```

```mathematica
In[1]:= WeightedAdjacencyMatrix[CycleGraph[3]]  (* unweighted: every present edge is 1 *)
```

### Notes

`WeightedAdjacencyMatrix[g]` is the dense adjacency matrix with each present
edge's `EdgeWeight` in place of a `1`; absent edges are `0`. An undirected edge
fills both symmetric cells, a directed edge only `M[tail][head]`.

For a graph with no weights every edge defaults to weight `1`, so the result
coincides with `AdjacencyMatrix[g]`. The row and column order is `VertexList`
order.
