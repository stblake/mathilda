### Worked examples

```mathematica
In[1]:= GraphDistanceMatrix[CycleGraph[4]]  (* entry (i, j) is the shortest-path length from i to j *)
```

```mathematica
In[1]:= GraphDistanceMatrix[PathGraph[{1, 2, 3, 4}], 1]  (* the second argument is a cutoff; longer distances become Infinity *)
```

```mathematica
In[1]:= GraphDistanceMatrix[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[2, 3]}]]  (* directed edges give an asymmetric matrix *)
```

```mathematica
In[1]:= GraphDistanceMatrix[Graph[{1, 2, 3}, {UndirectedEdge[1, 2]}]]  (* unreachable pairs are Infinity *)
```

```mathematica
In[1]:= Max[GraphDistanceMatrix[GridGraph[{3, 4}]]]  (* the largest entry is the diameter *)
```

### Notes

Rows are sources and columns are targets, in `VertexList` order. The diagonal is `0`. An unweighted graph gives exact integers; a graph with `EdgeWeight` gives machine reals.

When every entry is finite the result is a packed integer (or real) matrix. If any pair is unreachable or cut off by the distance bound, the matrix is returned as an ordinary list of lists with `Infinity` in those slots.
