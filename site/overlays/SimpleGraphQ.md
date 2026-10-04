### Worked examples

```mathematica
In[1]:= SimpleGraphQ[CompleteGraph[4]]  (* every valid Mathilda graph is simple *)
```

```mathematica
In[1]:= SimpleGraphQ[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]  (* a directed triangle is still simple *)
```

```mathematica
In[1]:= SimpleGraphQ[CycleGraph[6]]  (* an undirected cycle *)
```

```mathematica
In[1]:= SimpleGraphQ[5]  (* a non-graph argument is False *)
```

### Notes

A simple graph has no self-loops and no parallel (duplicate) edges. Mathilda's `Graph[...]`
constructor rejects both at build time, so every canonical graph is simple — `SimpleGraphQ`
is therefore `True` for any valid graph and reduces to a validity check. It coincides with
`LoopFreeGraphQ` for the same reason.

The predicate is memoized through the validated-graph cache, so it is `O(1)` after the
first query. A non-graph argument gives `False` rather than staying unevaluated.
