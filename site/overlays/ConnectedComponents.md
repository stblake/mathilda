### Worked examples

```mathematica
In[1]:= ConnectedComponents[Graph[{1 <-> 2, 3 <-> 4, 4 <-> 5}]]  (* undirected: largest component first *)
```

```mathematica
In[1]:= ConnectedComponents[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}]]  (* directed: strongly connected components *)
```

### Notes

For an undirected graph these are the connected components, listed largest first (ties by first
appearance). For a graph with directed edges they are the strongly connected components, listed
in reverse-topological order of the condensation — no edge runs from a component to a later one,
so sinks come first. Vertices inside a component keep `VertexList` order.

`ConnectedComponents[g, {v1, ...}]` keeps only the components containing one of the listed
vertices. `WeaklyConnectedComponents` ignores edge directions, and `StronglyConnectedComponents`
forces the directed interpretation even on an undirected graph.
