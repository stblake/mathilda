### Worked examples

```mathematica
In[1]:= DirectedGraphQ[Graph[{1, 2}, {1 -> 2}]]  (* every edge is a DirectedEdge *)
```

```mathematica
In[1]:= DirectedGraphQ[CycleGraph[3]]  (* an undirected graph is not directed *)
```

```mathematica
In[1]:= DirectedGraphQ[Graph[{1, 2}, {1 <-> 2}]]  (* an undirected edge, so False *)
```

### Notes

`DirectedGraphQ` is `True` only when the graph has at least one edge and all of
its edges are directed. An edgeless graph is treated as undirected, so it gives
`False` — and `DirectedGraphQ` and `UndirectedGraphQ` are never both `True` on
the same graph. A mixed graph (some directed, some undirected edges) gives
`False`; test for that with `MixedGraphQ`.
