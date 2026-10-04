### Worked examples

```mathematica
In[1]:= VertexCount[GraphDisjointUnion[CycleGraph[3], CycleGraph[4]]]  (* vertices never merge: 3 + 4 *)
```

```mathematica
In[1]:= EdgeCount[GraphDisjointUnion[CompleteGraph[3], CompleteGraph[3]]]  (* 3 + 3 edges, kept apart *)
```

```mathematica
In[1]:= VertexList[GraphDisjointUnion[CycleGraph[3], PathGraph[{1, 2}]]]  (* relabelled 1..n *)
```

### Notes

`GraphDisjointUnion` lays the graphs out on disjoint vertex sets even when they
share vertex names, relabelling everything to `1, ..., n` with the first graph's
vertices first. The result's vertex and edge counts are the plain sums of the
inputs'.

Contrast `GraphUnion`, which identifies equal vertices across its arguments.
Edge weights are dropped.
