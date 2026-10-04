### Worked examples

```mathematica
In[1]:= t = TuranGraph[7, 3];
```

```mathematica
In[1]:= {VertexCount[t], EdgeCount[t]}  (* parts of sizes 3, 2, 2 *)
```

```mathematica
In[1]:= EdgeCount[TuranGraph[6, 3]]  (* three parts of 2 is complete tripartite K(2,2,2) *)
```

### Notes

`TuranGraph[n, k]` is the complete `k`-partite graph on `n` vertices with parts as equal as
possible (the first `n mod k` parts get one extra vertex). By Turán's theorem it maximises the
edge count among `n`-vertex graphs that contain no `(k+1)`-clique.

The vertices are the integers `1..n` and the graph is always undirected. `TuranGraph[n, 2]` is
the balanced complete bipartite graph, and `TuranGraph[n, n]` is the complete graph.
