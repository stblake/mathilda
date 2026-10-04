### Worked examples

```mathematica
In[1]:= FindMinimumCut[PathGraph[{1, 2, 3}]]  (* one edge separates the path *)
```

```mathematica
In[1]:= FindMinimumCut[CycleGraph[4]]  (* a cycle must lose two edges to split *)
```

```mathematica
In[1]:= FindMinimumCut[CompleteGraph[4]]  (* isolating one vertex of K4 cuts three edges *)
```

```mathematica
In[1]:= FindMinimumCut[Graph[{1 -> 2, 2 -> 3, 3 -> 1, 1 -> 3}]]  (* directed: the source side is listed first *)
```

### Notes

`FindMinimumCut[g]` returns `{value, {part1, part2}}`, a global minimum edge cut:
the least total capacity of edges whose removal disconnects `g`, with the two
shores of that cut. Capacities come from `EdgeWeight` when `g` carries one, else
every edge has capacity `1`, so on an unweighted graph the value is the
edge-connectivity.

Undirected graphs use Nagamochi-Ibaraki (the Stoer-Wagner answer, computed by
edge contraction); graphs with a directed edge take the minimum over all vertex
pairs of the two directed max flows. The source side is listed first for a
directed graph, the side without the first vertex first for an undirected one.
Integer capacities give an exact integer value, rational or real ones a real; a
graph on fewer than two vertices, or a negative/symbolic capacity, leaves the
call unevaluated.
