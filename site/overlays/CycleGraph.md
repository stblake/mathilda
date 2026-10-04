### Worked examples

```mathematica
In[1]:= EdgeCount[CycleGraph[5]]  (* n edges when n >= 3 *)
```

```mathematica
In[1]:= EdgeList[CycleGraph[4]]  (* the rim, closing with 4 <-> 1 *)
```

```mathematica
In[1]:= VertexDegree[CycleGraph[5]]  (* every vertex of a cycle has degree two *)
```

```mathematica
In[1]:= VertexList[CycleGraph[6]]  (* the integer vertices 1..n *)
```

```mathematica
In[1]:= EdgeCount[CycleGraph[2]]  (* n <= 2 drops the duplicate wrap edge *)
```

### Notes

`CycleGraph[n]` is the undirected cycle `C_n` on the integer vertices `1..n`,
closing the path `1-2-...-n` with the edge `n <-> 1`. For `n >= 3` it has exactly
`n` edges and every vertex has degree two; `n = 2` is a single edge (the wrap edge
would duplicate it) and `n <= 1` is edgeless. The result is an opaque `Graph`, read
through accessors such as `EdgeList`, `VertexDegree`, and `EdgeCount`.
