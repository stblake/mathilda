### Worked examples

```mathematica
In[1]:= VertexList[VertexAdd[PathGraph[{1, 2, 3}], 7]]  (* a single new vertex goes at the end *)
```

```mathematica
In[1]:= VertexList[VertexAdd[PathGraph[{1, 2, 3}], {4, 5, 2, 4}]]  (* existing and repeated vertices are ignored *)
```

```mathematica
In[1]:= EdgeCount[VertexAdd[CycleGraph[4], {a, b}]]  (* new vertices are isolated, the edges are untouched *)
```

```mathematica
In[1]:= VertexCount[VertexAdd[Graph[{1 -> 2}], {x, y, z}]]  (* works on directed graphs too *)
```

### Notes

Added vertices are appended after the existing ones, in the order given, and never carry edges. A vertex already present is skipped, so `VertexAdd` is idempotent.
