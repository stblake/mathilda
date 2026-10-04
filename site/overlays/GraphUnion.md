### Worked examples

```mathematica
In[1]:= VertexCount[GraphUnion[CycleGraph[3], PathGraph[{3, 4, 5}]]]  (* the vertices are the set-union 1..5 *)
```

```mathematica
In[1]:= EdgeList[GraphUnion[CycleGraph[3], PathGraph[{3, 4, 5}]]]  (* every distinct edge of either graph *)
```

```mathematica
In[1]:= EdgeList[GraphUnion[PathGraph[{1, 2, 3}], PathGraph[{2, 3, 4}]]]  (* the shared edge 2 <-> 3 appears once *)
```

```mathematica
In[1]:= EdgeCount[GraphUnion[CompleteGraph[3], CompleteGraph[3]]]  (* a graph unioned with itself is itself *)
```

### Notes

`GraphUnion[g1, g2, ...]` is the graph whose vertices are the union of all the
graphs' vertices (in canonical `Sort` order) and whose edges are every **distinct**
edge of any of them — an undirected edge and its reversal count as one, so a shared
edge appears a single time. `GraphUnion[g]` is just `g`. The result is an opaque
`Graph`, read through `VertexCount`, `EdgeCount`, `EdgeList`, and the other
accessors.

Edge order matches Mathematica 15 (first-appearance for a uniform direction,
canonical order when directed and undirected edges are mixed), and edge **weights
are dropped** — the union is an unweighted simple graph.
