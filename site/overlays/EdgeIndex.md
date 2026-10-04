### Worked examples

```mathematica
In[1]:= EdgeIndex[CycleGraph[4], 2 <-> 3]  (* position of this edge in EdgeList *)
```

```mathematica
In[1]:= EdgeIndex[CompleteGraph[3], 1 <-> 3]  (* K3 lists {1<->2, 1<->3, 2<->3} *)
```

```mathematica
In[1]:= EdgeIndex[CycleGraph[5], {1 <-> 2, 3 <-> 4}]  (* a list of edges gives a list of positions *)
```

### Notes

`EdgeIndex[g, e]` is the inverse of indexing `EdgeList[g]`: it reports where `e`
sits in the edge list, counting from 1. An undirected edge matches in either
orientation, and the sugar `u -> v` / `u <-> v` is accepted alongside explicit
`DirectedEdge` / `UndirectedEdge`.

An edge that is not in the graph leaves the expression unevaluated rather than
returning `0` or a `Missing` wrapper.
