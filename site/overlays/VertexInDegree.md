### Worked examples

```mathematica
In[1]:= VertexInDegree[Graph[{1 -> 2, 1 -> 3, 2 -> 3}]]  (* in-degrees in vertex order *)
```

```mathematica
In[1]:= VertexInDegree[Graph[{1 -> 2, 1 -> 3, 2 -> 3}], 3]  (* the in-degree of a single vertex *)
```

```mathematica
In[1]:= VertexInDegree[CycleGraph[4]]  (* undirected: in-degree equals total degree *)
```

### Notes

A directed edge `a -> b` contributes to the in-degree of `b` only. An undirected
edge is incident to both endpoints, so for a purely undirected graph in-degree
equals out-degree equals total degree.

`VertexInDegree[g]` gives the list of in-degrees in canonical vertex order;
`VertexInDegree[g, v]` gives the single in-degree of `v`, and leaves itself
unevaluated when `v` is not a vertex. A `Hypergraph` has only a total degree, so
the in-degree form does not apply to one.
