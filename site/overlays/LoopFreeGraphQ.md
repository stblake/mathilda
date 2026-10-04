### Worked examples

```mathematica
In[1]:= LoopFreeGraphQ[CycleGraph[3]]  (* no edge joins a vertex to itself *)
```

```mathematica
In[1]:= LoopFreeGraphQ[CompleteGraph[4]]  (* complete graphs have no self-loops *)
```

```mathematica
In[1]:= LoopFreeGraphQ[Graph[{1 -> 2, 2 -> 3}]]  (* a directed path is loop-free *)
```

```mathematica
In[1]:= LoopFreeGraphQ["x"]  (* a non-graph argument is False *)
```

### Notes

A loop-free graph has no self-loops — no edge from a vertex to itself. Mathilda's `Graph[...]`
constructor refuses self-loops, so every valid graph is loop-free and the predicate reduces
to a validity check. For the same structural reason it agrees with `SimpleGraphQ` on every
input.

The answer is memoized through the validated-graph cache (`O(1)` after the first query),
and a non-graph argument yields `False` rather than remaining unevaluated.
