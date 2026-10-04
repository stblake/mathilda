### Worked examples

```mathematica
In[1]:= FindPath[CycleGraph[6], 1, 4]  (* first path found by depth-first search *)
```

```mathematica
In[1]:= FindPath[CycleGraph[6], 1, 4, {3}, All]  (* all paths of exactly three edges *)
```

```mathematica
In[1]:= FindPath[CompleteGraph[4], 1, 4, Infinity, All]  (* every simple path, shortest first *)
```

```mathematica
In[1]:= FindPath[Graph[{1 -> 2, 2 -> 3}], 3, 1]  (* no directed path against the arrows *)
```

```mathematica
In[1]:= FindPath[Graph[{1 <-> 2, 3 <-> 4}], 1, 4]  (* disconnected components give no path *)
```

```mathematica
In[1]:= FindPath[CycleGraph[6], 1, 1]  (* start equal to end gives an empty result *)
```

### Notes

A path is a list of vertices and the result is a list of paths, `{}` if none exists. Length counts edges, with the same `kspec` forms as `FindCycle`, and a fifth argument bounds the number of paths (or `All`); multiple paths are listed shortest first.

The single-path search is a linear-time depth-first search with neighbours taken in `EdgeList` order. Enumeration is bounded by a step budget, past which the call stays unevaluated.
