### Worked examples

```mathematica
In[1]:= EdgeList[StarGraph[5]]  (* hub 1 joined to four leaves *)
```

```mathematica
In[1]:= VertexDegree[StarGraph[6]]  (* the hub has degree n - 1, each leaf degree 1 *)
```

```mathematica
In[1]:= EdgeCount[StarGraph[1]]  (* a single vertex has no edges *)
```

```mathematica
In[1]:= AdjacencyMatrix[StarGraph[4]]  (* only the first row and column carry ones *)
```

```mathematica
In[1]:= GraphQ[StarGraph[8]]  (* the result is an ordinary Graph expression *)
```

```mathematica
In[1]:= VertexCount[StarGraph[100]]  (* cheap even for large n, with n - 1 edges *)
```

### Notes

`StarGraph[n]` has `n` vertices with vertex 1 as the hub joined to each of `2..n`; it is the complete bipartite graph `K(1, n-1)`. `StarGraph[1]` is one isolated vertex.

The argument must be a non-negative integer; anything else leaves the call unevaluated.
