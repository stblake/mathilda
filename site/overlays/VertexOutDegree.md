### Worked examples

```mathematica
In[1]:= VertexOutDegree[Graph[{1, 2, 3}, {1 -> 2, 1 -> 3}]]  (* vertex 1 has two out-edges *)
```

```mathematica
In[1]:= VertexOutDegree[Graph[{1, 2, 3}, {1 -> 2, 1 -> 3}], 1]  (* the out-degree of a single vertex *)
```

```mathematica
In[1]:= VertexOutDegree[CycleGraph[4]]  (* an undirected edge counts for out-degree too *)
```

### Notes

`VertexOutDegree[g]` gives the number of edges leaving each vertex, in
`VertexList` order; `VertexOutDegree[g, v]` the value for one vertex. A directed
edge contributes to the out-degree of its tail only; an undirected edge is
incident to both ends, so it contributes to the out-degree of each — which is
why on an undirected graph the out-degree equals the ordinary degree.

The companions are `VertexInDegree` and `VertexDegree` (the total).
