### Worked examples

```mathematica
In[1]:= PlanarGraphQ[CompleteGraph[4]]  (* k4 draws in the plane without crossings *)
```

```mathematica
In[1]:= PlanarGraphQ[CompleteGraph[5]]  (* k5 is the first forbidden minor *)
```

```mathematica
In[1]:= PlanarGraphQ[CompleteGraph[{3, 3}]]  (* so is the utility graph k3,3 *)
```

```mathematica
In[1]:= PlanarGraphQ[GridGraph[{3, 3}]]  (* grids are planar *)
```

```mathematica
In[1]:= PlanarGraphQ[PetersenGraph[]]  (* ten vertices, fifteen edges, but not planar *)
```

```mathematica
In[1]:= PlanarGraphQ[Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}]]  (* direction is ignored *)
```

```mathematica
In[1]:= PlanarGraphQ[x]  (* a non-graph gives false *)
```

### Notes

The test is the linear-time left-right criterion, so it is cheap even on very large sparse graphs. It decides planarity only; it does not return an embedding or a Kuratowski subgraph.

Edge direction, anti-parallel pairs and loops are ignored: the test runs on the underlying simple undirected graph. Graphs with more than `3n - 6` edges are rejected immediately by the Euler bound.
