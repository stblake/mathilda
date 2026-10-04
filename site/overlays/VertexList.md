### Worked examples

```mathematica
In[1]:= VertexList[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]  (* vertices derived from the edges *)
```

```mathematica
In[1]:= VertexList[PathGraph[{a, b, c}]]  (* vertices may be any expressions *)
```

### Notes

`VertexList[g]` returns the vertices of `g` in canonical order. When a graph is built from edges
alone, the vertices are derived from the edge endpoints in order of first appearance.

Vertices are arbitrary expressions, not just integers; a path on symbols `a`, `b`, `c` lists
them unchanged.
