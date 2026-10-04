### Worked examples

```mathematica
In[1]:= EdgeWeightedGraphQ[Graph[{1 <-> 2}, EdgeWeight -> {2}]]  (* one weighted undirected edge *)
```

```mathematica
In[1]:= EdgeWeightedGraphQ[Graph[{1 -> 2, 2 -> 3}, EdgeWeight -> {3, 4}]]  (* a weighted directed path *)
```

```mathematica
In[1]:= EdgeWeightedGraphQ[CycleGraph[4]]  (* no EdgeWeight, so False *)
```

```mathematica
In[1]:= EdgeWeightedGraphQ[0]  (* a non-graph argument is False *)
```

### Notes

`EdgeWeightedGraphQ[g]` is `True` exactly when `g` carries an `EdgeWeight` option — the
same test as `WeightedGraphQ`, since Mathilda associates weights only with edges. It checks
for the presence of the weight list, not that the weights are numeric.

A graph built without `EdgeWeight` is unweighted (each edge's weight defaults to `1` in the
accessors) and gives `False`, as does any non-graph argument.
