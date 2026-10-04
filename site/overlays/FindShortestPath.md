### Worked examples

```mathematica
In[1]:= FindShortestPath[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4}], 1, 4]  (* follows edge directions *)
```

```mathematica
In[1]:= FindShortestPath[Graph[{1 -> 2, 3 -> 4}], 1, 4]  (* no path exists *)
```

### Notes

The result is the path from `s` to `t` as a list of vertices; it is `{}` when no path exists.
On a directed graph the path respects edge directions.

A graph that carries non-negative numeric `EdgeWeight`s is traversed with Dijkstra's algorithm,
so the path minimises total weight; an unweighted graph (or one with symbolic or negative
weights) falls back to a breadth-first search that minimises the hop count. `GraphDistance`
gives the corresponding path length.
