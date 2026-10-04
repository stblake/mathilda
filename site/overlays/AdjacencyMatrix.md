### Worked examples

```mathematica
In[1]:= AdjacencyMatrix[CycleGraph[4]]  (* symmetric for an undirected graph *)
```

```mathematica
In[1]:= AdjacencyMatrix[Graph[{1 -> 2, 2 -> 3}]]  (* a directed edge sets only M[a][b] *)
```

```mathematica
In[1]:= AdjacencyMatrix[StarGraph[4]]  (* hub row and column *)
```

```mathematica
In[1]:= Tr[MatrixPower[AdjacencyMatrix[CompleteGraph[3]], 3]]  (* six times the number of triangles *)
```

```mathematica
In[1]:= Det[AdjacencyMatrix[CycleGraph[4]]]  (* feeds linear algebra directly *)
```

```mathematica
In[1]:= AdjacencyMatrix[Graph[{a, b, c}, {a <-> c}]]  (* rows follow the vertex list order *)
```

### Notes

`AdjacencyMatrix` returns a dense `List` of `List`s of `0` and `1`, with rows in `VertexList` order. An undirected edge fills both `M[a][b]` and `M[b][a]`, so such graphs give a symmetric matrix. A directed edge fills only `M[a][b]`. Because the result is an ordinary matrix, `Det`, `Tr`, `MatrixPower` and `Eigenvalues` apply directly.
