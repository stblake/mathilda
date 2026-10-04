### Worked examples

```mathematica
In[1]:= FindMaximumFlow[Graph[{1, 2, 3}, {1 -> 2, 2 -> 3, 1 -> 3}, EdgeCapacity -> {2, 1, 3}], 1, 3]  (* direct 3 plus 1 through vertex 2 *)
```

```mathematica
In[1]:= FindMaximumFlow[CompleteGraph[4], 1, 2]  (* unit capacities: three edge-disjoint paths *)
```

```mathematica
In[1]:= FindMaximumFlow[Graph[{1, 2, 3, 4}, {1 -> 2, 1 -> 3, 2 -> 4, 3 -> 4}, EdgeCapacity -> {3, 2, 2, 3}], 1, 4]  (* bottlenecked by the middle edges *)
```

### Notes

The flow value is capacity-limited: capacities come from the `EdgeCapacity`
option (`EdgeList` order), defaulting to `1` per edge when none is given, and
`EdgeWeight` plays no part. With integer capacities the answer is exact.

A fourth argument requests `"FlowValue"`, `"FlowMatrix"`, or `"EdgeList"`.
Sources and sinks may each be a list, and `VertexCapacity` bounds the flow
through individual vertices. An undirected edge can carry flow in either
direction; a directed edge only forwards.
