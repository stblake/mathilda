### Worked examples

```mathematica
In[1]:= KirchhoffMatrix[CycleGraph[4]]  (* degree 2 on the diagonal, -1 for each neighbour *)
```

```mathematica
In[1]:= KirchhoffMatrix[StarGraph[4]]  (* the hub has degree 3, the leaves degree 1 *)
```

```mathematica
In[1]:= KirchhoffMatrix[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[2, 3]}]]  (* a directed edge fills only one off-diagonal entry *)
```

```mathematica
In[1]:= Total[KirchhoffMatrix[CompleteGraph[5]], 2]  (* undirected rows sum to zero *)
```

```mathematica
In[1]:= Det[KirchhoffMatrix[CycleGraph[4]]]  (* the Laplacian is singular *)
```

### Notes

The result is `D - A`, where `D` holds the number of edges incident to each vertex and `A` is the adjacency matrix. Edge weights are ignored.

Mathematica returns a `SparseArray`; Mathilda returns the equivalent dense integer matrix, which is what `Normal` of Mathematica's answer gives.
