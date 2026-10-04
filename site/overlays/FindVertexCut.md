### Worked examples

```mathematica
In[1]:= FindVertexCut[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 2 <-> 4, 4 <-> 5}]]  (* one cut vertex suffices *)
```

```mathematica
In[1]:= FindVertexCut[CompleteGraph[4]]  (* a complete graph needs n-1 vertices removed *)
```

### Notes

`FindVertexCut[g]` gives a smallest vertex set whose removal disconnects the underlying
undirected graph; its length equals `VertexConnectivity[g]`. For a complete graph on `n`
vertices this is `n-1` vertices (you cannot disconnect a clique without deleting all but one).

`FindVertexCut[g, s, t]` gives a smallest set separating `s` from `t`, and returns `{}` when `s`
and `t` are adjacent.
