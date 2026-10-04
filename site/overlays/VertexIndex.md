### Worked examples

```mathematica
In[1]:= VertexIndex[CycleGraph[5], 3]  (* position of a vertex in VertexList *)
```

```mathematica
In[1]:= VertexIndex[Graph[{a, b, c}, {a <-> b, b <-> c}], c]  (* symbolic vertices work too *)
```

```mathematica
In[1]:= VertexIndex[Graph[{a, b, c}, {a <-> b, b <-> c}], {c, a}]  (* a list gives a list of positions *)
```

```mathematica
In[1]:= VertexIndex[CycleGraph[5], 5]  (* the last vertex of a pentagon *)
```

### Notes

The index is 1-based and follows `VertexList`, which need not be sorted. Every vertex in a list argument must exist, otherwise the whole call is left unevaluated.
