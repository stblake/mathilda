### Worked examples

```mathematica
In[1]:= FindEdgeCut[PathGraph[{1, 2, 3}]]  (* a single bridge is enough *)
```

```mathematica
In[1]:= FindEdgeCut[CycleGraph[4]]  (* two edges, here both incident to one vertex *)
```

```mathematica
In[1]:= FindEdgeCut[CompleteGraph[4], 1, 2]  (* a minimum 1-2 edge cut *)
```

### Notes

`FindEdgeCut[g]` returns the edges of a global minimum cut; `FindEdgeCut[g, s,
t]` the edges of a minimum cut separating `s` from `t`. The number of edges (or
their total `EdgeWeight`) equals `EdgeConnectivity` of the same arguments.

Edges are returned in `EdgeList` order. The `s-t` cut is the one closest to
`s`; adjacent `s` and `t` on an unweighted graph can only be separated by
removing every edge between them.
