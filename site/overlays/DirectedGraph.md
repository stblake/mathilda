### Worked examples

```mathematica
In[1]:= EdgeList[DirectedGraph[Graph[{1 <-> 2, 2 <-> 3}]]]  (* each undirected edge becomes a pair of arcs *)
```

```mathematica
In[1]:= EdgeList[DirectedGraph[Graph[{1 <-> 2, 2 <-> 3}], "Acyclic"]]  (* oriented low-to-high: a DAG *)
```

### Notes

By default `DirectedGraph[g]` replaces each undirected edge `u <-> v` by the two arcs `u -> v`
and `v -> u` (weights duplicated), so the result has the same reachability as `g`.

`DirectedGraph[g, "Acyclic"]` instead orients each undirected edge from the earlier to the later
vertex in `VertexList`, producing a DAG. A graph that is already fully directed is returned
unchanged. The result is a canonical `Graph` — inspect it with `EdgeList`.
