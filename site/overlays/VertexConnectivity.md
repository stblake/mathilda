### Worked examples

```mathematica
In[1]:= VertexConnectivity[CompleteGraph[4]]  (* a complete graph on n needs n-1 removed *)
```

```mathematica
In[1]:= VertexConnectivity[PathGraph[5]]  (* removing one interior vertex breaks a path *)
```

```mathematica
In[1]:= VertexConnectivity[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 2 <-> 4, 4 <-> 5}]]  (* vertex 4 is a cut vertex *)
```

### Notes

The result is the smallest number of vertices whose removal disconnects the graph (or reduces it
to a single vertex). A complete graph on `n` vertices has connectivity `n-1`; a graph that is
already disconnected, or has at most one vertex, has connectivity `0`.

`FindVertexCut` returns an actual minimum separating set of that size.
