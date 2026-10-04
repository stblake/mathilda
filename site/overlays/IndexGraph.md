### Worked examples

```mathematica
In[1]:= VertexList[IndexGraph[Graph[{a, b, c}, {a <-> b, b <-> c}]]]  (* vertices become 1, 2, 3 *)
```

```mathematica
In[1]:= EdgeList[IndexGraph[Graph[{a, b, c}, {a <-> b, b <-> c}], 0]]  (* an offset of 0 starts at 0 *)
```

```mathematica
In[1]:= EdgeList[IndexGraph[Graph[{x -> y, y -> z}]]]  (* directed edges keep their direction *)
```

```mathematica
In[1]:= EdgeList[IndexGraph[CycleGraph[4], 10]]  (* offset labels r, r+1, ... *)
```

### Notes

Vertex number `i` in `VertexList` is renamed `r + i - 1`, with `r = 1` by default, so the structure of the graph is unchanged and only the labels differ. The offset must be a machine integer.
