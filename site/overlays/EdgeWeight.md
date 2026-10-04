### Worked examples

```mathematica
In[1]:= EdgeWeight[CycleGraph[3]]  (* an unweighted graph defaults every weight to 1 *)
```

```mathematica
In[1]:= EdgeWeight[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}, EdgeWeight -> {5, 7}]]  (* the stored weights, in EdgeList order *)
```

### Notes

`EdgeWeight[g]` always returns a weight for every edge, in `EdgeList` order. A
graph built without weights reports `1` for each edge, matching the Wolfram
convention that an unweighted edge has weight `1` — so `WeightedAdjacencyMatrix`
and the weighted shortest-path routines have a defined weight to use on any
graph.

The weight list is matched to the edges by position: `EdgeWeight[g][[k]]` is the
weight of `EdgeList[g][[k]]`.
