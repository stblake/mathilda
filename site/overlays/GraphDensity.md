### Worked examples

```mathematica
In[1]:= GraphDensity[CompleteGraph[4]]  (* a complete graph has density 1 *)
```

```mathematica
In[1]:= GraphDensity[CycleGraph[5]]  (* five edges out of ten possible *)
```

```mathematica
In[1]:= GraphDensity[Graph[{1 -> 2, 2 -> 3}]]  (* directed arcs each count once *)
```

```mathematica
In[1]:= GraphDensity[CompleteGraph[{2, 3}]]  (* a complete bipartite graph *)
```

### Notes

The density is the number of edges present divided by the number that could be
present: `(directed arcs + 2 × undirected edges) / (n(n−1))`. An undirected edge
counts twice because it fills both ordered endpoint pairs, so a complete undirected
graph scores exactly `1` and a sparse graph scores near `0`.

The result is an exact rational, computed from the vertex and edge counts alone
with no traversal. A graph on fewer than two vertices is left unevaluated, since
`n(n−1)` would be zero and the density is undefined.
