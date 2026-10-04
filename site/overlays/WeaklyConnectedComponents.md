### Worked examples

```mathematica
In[1]:= WeaklyConnectedComponents[Graph[{1 -> 2, 3 -> 4}]]  (* directions ignored: two components *)
```

### Notes

The weakly connected components are the components of the underlying undirected graph: two
vertices are together when they are joined by a path that ignores edge directions. They are
listed largest first (ties by first appearance), with vertices inside a component in
`VertexList` order.

`WeaklyConnectedComponents[g, {v1, ...}]` keeps only the components containing one of the listed
vertices. For a directed graph with a cycle, the weak components can be strictly coarser than
the strongly connected ones.
