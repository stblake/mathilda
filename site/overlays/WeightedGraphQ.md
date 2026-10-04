### Worked examples

```mathematica
In[1]:= WeightedGraphQ[Graph[{1 -> 2, 2 -> 3}, EdgeWeight -> {5, 7}]]  (* carries an EdgeWeight list *)
```

```mathematica
In[1]:= WeightedGraphQ[Graph[{1 <-> 2}, EdgeWeight -> {2}]]  (* an undirected weighted edge *)
```

```mathematica
In[1]:= WeightedGraphQ[CompleteGraph[3]]  (* a generator graph carries no weights *)
```

```mathematica
In[1]:= WeightedGraphQ[CycleGraph[4]]  (* unweighted, so False *)
```

### Notes

`WeightedGraphQ[g]` is `True` exactly when `g` is a valid graph built with an `EdgeWeight`
option. It coincides with `EdgeWeightedGraphQ` — Mathilda has no vertex weights, so a graph
is weighted iff its edges are. The test checks only for the *presence* of the weight list,
not that the weights are numeric or usable by a shortest-path algorithm.

An unweighted graph (every accessor then treats each edge's weight as `1`) gives `False`,
as does any non-graph argument.
