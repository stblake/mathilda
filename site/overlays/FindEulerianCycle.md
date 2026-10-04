### Worked examples

```mathematica
In[1]:= FindEulerianCycle[CycleGraph[4]]  (* every vertex has even degree, so a closed Euler tour exists *)
```

```mathematica
In[1]:= FindEulerianCycle[CompleteGraph[3]]  (* the triangle is Eulerian *)
```

```mathematica
In[1]:= FindEulerianCycle[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]  (* directed: in-degree equals out-degree at each vertex *)
```

```mathematica
In[1]:= FindEulerianCycle[PathGraph[{1, 2, 3}]]  (* odd-degree endpoints, so no Euler cycle -- an empty result *)
```

```mathematica
In[1]:= FindEulerianCycle[CompleteGraph[4]]  (* every vertex has odd degree 3, so again none exists *)
```

### Notes

`FindEulerianCycle[g]` returns `{cycle}`, where the cycle is the list of edges of a
closed walk that uses every edge of `g` exactly once, or `{}` when no such walk
exists. A graph is Eulerian when it is connected and every vertex has even degree
(undirected) or equal in- and out-degree (directed); the head checks that first and
returns `{}` immediately otherwise.

The tour is built by Hierholzer's algorithm, starting from the first vertex (in
`VertexList` order) that has an incident edge and taking edges in `EdgeList` order,
so the walk is deterministic. Only `FindEulerianCycle[g]` and `FindEulerianCycle[g,
1]` are supported; a mixed directed/undirected graph, or a request for more than
one tour, is left unevaluated.
