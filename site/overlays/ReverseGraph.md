### Worked examples

```mathematica
In[1]:= EdgeList[ReverseGraph[Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}]]]  (* each arc flips direction *)
```

```mathematica
In[1]:= EdgeList[ReverseGraph[CycleGraph[3]]]  (* undirected edges are unchanged *)
```

### Notes

`ReverseGraph[g]` reverses the direction of every directed edge; undirected
edges, the edge order, and any `EdgeWeight` are kept. An undirected graph is
returned unchanged.

Reversing twice recovers the original graph.
