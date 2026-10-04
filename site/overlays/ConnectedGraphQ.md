### Worked examples

```mathematica
In[1]:= ConnectedGraphQ[PathGraph[{1, 2, 3}]]  (* a path is connected *)
```

```mathematica
In[1]:= ConnectedGraphQ[Graph[{1, 2, 3}, {1 <-> 2}]]  (* vertex 3 is isolated *)
```

```mathematica
In[1]:= ConnectedGraphQ[CycleGraph[4]]  (* a cycle is connected *)
```

### Notes

For an undirected graph this is ordinary connectivity: one component containing
every vertex. For a graph with directed edges it is *strong* connectivity —
every vertex must reach every other along the edge directions — matching
Mathematica, where a directed graph is "connected" only when strongly connected.

The empty graph (no vertices) is not connected.
