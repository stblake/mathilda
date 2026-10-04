### Worked examples

```mathematica
In[1]:= IncidenceMatrix[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]  (* oriented: -1 tail, +1 head *)
```

```mathematica
In[1]:= MatrixForm[IncidenceMatrix[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]]
```

### Notes

The result is a `|V| x |E|` matrix: rows index vertices in `VertexList` order, columns index
edges in `EdgeList` order. A directed edge is oriented — `-1` in its tail row and `+1` in its
head row; an undirected edge puts `+1` in both endpoint rows.

The matrix is a plain nested list, so it feeds straight into `MatrixForm` or the linear algebra
builtins.
