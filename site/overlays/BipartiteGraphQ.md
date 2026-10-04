### Worked examples

```mathematica
In[1]:= BipartiteGraphQ[CycleGraph[4]]  (* an even cycle is bipartite *)
```

```mathematica
In[1]:= BipartiteGraphQ[CycleGraph[5]]  (* an odd cycle is not *)
```

```mathematica
In[1]:= BipartiteGraphQ[CompleteGraph[{2, 3}]]  (* every complete bipartite graph qualifies *)
```

```mathematica
In[1]:= BipartiteGraphQ[CompleteGraph[4]]  (* K4 contains a triangle, so False *)
```

### Notes

A graph is bipartite exactly when it has no odd cycle. The test is a BFS
2-colouring over each connected component; direction is ignored, so a directed
graph is tested as its underlying undirected graph.

The answer is cached on the graph object, so repeated bipartiteness queries on
the same `Graph` cost nothing after the first. A non-graph argument gives
`False` rather than staying unevaluated.
