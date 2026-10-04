### Worked examples

```mathematica
In[1]:= VertexCount[CompleteGraph[6]]  (* K6 has six vertices *)
```

```mathematica
In[1]:= VertexCount[CycleGraph[4]]
```

```mathematica
In[1]:= VertexCount[Graph[{a, b, c}, {a <-> b}]]  (* vertices may be arbitrary symbols *)
```

### Notes

`VertexCount[g]` is the number of vertices, read directly off the canonical
graph form. Vertices can be any expressions, not just integers, and an isolated
vertex (one with no incident edge) still counts.

The parallel reader for edges is `EdgeCount`.
