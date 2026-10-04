### Worked examples

```mathematica
In[1]:= VertexDegree[StarGraph[5]]  (* one degree per vertex, in vertex order *)
```

```mathematica
In[1]:= VertexDegree[CycleGraph[6], 1]  (* the degree of a single vertex *)
```

```mathematica
In[1]:= VertexDegree[CompleteGraph[4]]  (* every vertex meets the other three *)
```

```mathematica
In[1]:= VertexDegree[Graph[{1 -> 2, 1 -> 3, 3 -> 1}]]  (* directed edges count at both ends *)
```

```mathematica
In[1]:= Total[VertexDegree[PetersenGraph[]]] == 2 EdgeCount[PetersenGraph[]]  (* the handshake lemma *)
```

```mathematica
In[1]:= VertexDegree[Graph[{a, b, c}, {a <-> b}]]  (* an isolated vertex has degree zero *)
```

### Notes

`VertexDegree[g]` lists degrees in `VertexList` order, and `VertexDegree[g, v]` gives one vertex. In a directed graph an edge counts at both endpoints, so total degree is in-degree plus out-degree. `VertexInDegree` and `VertexOutDegree` share the same code. A `v` that is not a vertex leaves the call unevaluated.
