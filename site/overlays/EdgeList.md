### Worked examples

```mathematica
In[1]:= EdgeList[CompleteGraph[3]]  (* undirected edges, in canonical order *)
```

```mathematica
In[1]:= EdgeList[Graph[{1 -> 2, 2 -> 3}]]  (* directed edges keep their arrows *)
```

```mathematica
In[1]:= EdgeCount[CycleGraph[4]] == Length[EdgeList[CycleGraph[4]]]  (* the count is the list length *)
```

### Notes

`EdgeList[g]` returns the graph's edges in the canonical two-argument form:
`UndirectedEdge[u, v]` (printed `u <-> v`) and `DirectedEdge[u, v]` (printed
`u -> v`). Edge sugar given to the constructor is already normalised to these
heads, so `EdgeList` round-trips the edge half of `Graph[...]`.

A `Hypergraph` argument returns its hyperedges; any other non-graph argument
leaves the call unevaluated.
