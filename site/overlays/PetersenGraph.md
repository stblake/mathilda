### Worked examples

```mathematica
In[1]:= VertexCount[PetersenGraph[]]  (* the Petersen graph has 10 vertices *)
```

```mathematica
In[1]:= EdgeCount[PetersenGraph[]]  (* and 15 edges *)
```

```mathematica
In[1]:= VertexDegree[PetersenGraph[]]  (* it is cubic: every vertex has degree 3 *)
```

```mathematica
In[1]:= VertexList[PetersenGraph[3, 1]]  (* the generalized form GP(n, k) has 2n vertices labelled 1..2n *)
```

```mathematica
In[1]:= VertexCount[PetersenGraph[7, 2]]  (* a larger generalized Petersen graph *)
```

### Notes

`PetersenGraph[]` is the classic Petersen graph `GP(5, 2)`: an inner pentagram
joined by spokes to an outer pentagon, 3-regular on 10 vertices. `PetersenGraph[n,
k]` gives the generalized Petersen graph `GP(n, k)` on `2n` vertices, with inner
vertices `1 … n` joined at step `k` and outer vertices `n+1 … 2n` forming a cycle.

The construction requires that `k` not be a multiple of `n`, so the inner edges are
real chords rather than self-loops. The result is an opaque `Graph` object; query a
property (`VertexCount`, `EdgeCount`, `VertexDegree`, `VertexList`) to inspect it.
