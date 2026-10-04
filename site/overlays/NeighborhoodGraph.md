### Worked examples

```mathematica
In[1]:= VertexList[NeighborhoodGraph[CycleGraph[6], 1]]  (* vertex 1 and its two cycle neighbours *)
```

```mathematica
In[1]:= VertexList[NeighborhoodGraph[CycleGraph[6], 1, 2]]  (* the ball of radius two around vertex 1 *)
```

```mathematica
In[1]:= VertexList[NeighborhoodGraph[PathGraph[{1, 2, 3, 4, 5}], 3, 1]]  (* the centre first, then its ball in order *)
```

```mathematica
In[1]:= EdgeCount[NeighborhoodGraph[CompleteGraph[5], 1]]  (* in K5 one step reaches every vertex *)
```

### Notes

`NeighborhoodGraph[g, v, k]` is the subgraph of `g` induced by `v` and every vertex
within graph distance `k` of it (`k` defaults to `1`, and `Infinity` is allowed).
Edge direction is ignored when measuring distance. The centre may instead be a list
of vertices, and non-vertices among them are quietly skipped. The result is an
opaque `Graph`, read through `VertexList`, `EdgeList`, `EdgeCount`, and the others.

Vertices come out centres-first, then each centre's newly reached ball in
`VertexList` order — so `NeighborhoodGraph[CycleGraph[6], 1]` lists `1` before its
neighbours `2` and `6`. With `k` large enough the result is the whole connected
component of the centres.
