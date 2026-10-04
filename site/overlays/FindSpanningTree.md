### Worked examples

```mathematica
In[1]:= EdgeList[FindSpanningTree[CycleGraph[4]]]  (* drops one edge to break the cycle *)
```

```mathematica
In[1]:= EdgeCount[FindSpanningTree[CompleteGraph[5]]]  (* a spanning tree has n-1 edges *)
```

### Notes

The result is a `Graph` on the same vertices whose edges form a spanning tree
(or forest, if the graph is disconnected). A tree on `n` vertices has `n - 1`
edges.

When the graph carries `EdgeWeight`, the tree is a *minimum* spanning tree —
Kruskal for undirected graphs, Chu-Liu/Edmonds for directed ones — and weights
are compared exactly, so an integer and the nearest machine real are never
confused. A mixed graph is left unevaluated.
