### Worked examples

```mathematica
In[1]:= EdgeList[UndirectedGraph[Graph[{1, 2, 3}, {1 -> 2, 2 -> 3, 3 -> 1}]]]  (* a directed triangle becomes undirected *)
```

```mathematica
In[1]:= UndirectedGraphQ[UndirectedGraph[Graph[{1, 2}, {1 -> 2}]]]  (* direction dropped *)
```

```mathematica
In[1]:= EdgeList[UndirectedGraph[Graph[{1, 2}, {1 -> 2, 2 -> 1}]]]  (* the two arcs merge into one edge *)
```

### Notes

`UndirectedGraph[g]` drops direction from every edge. A pair of opposite arcs
`u -> v` and `v -> u` collapses into the single undirected edge `u <-> v`; when
the graph is weighted, the merged edge's weight is the sum of the two arcs'
weights.

Edges come out oriented and ordered by `VertexList` position (upper triangle,
row-major). An already-undirected graph is returned unchanged. The inverse-
direction operation is `DirectedGraph`.
