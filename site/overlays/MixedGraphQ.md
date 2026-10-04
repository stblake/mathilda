### Worked examples

```mathematica
In[1]:= MixedGraphQ[Graph[{1, 2, 3}, {1 -> 2, 2 <-> 3}]]  (* one directed and one undirected edge *)
```

```mathematica
In[1]:= MixedGraphQ[CycleGraph[3]]  (* all undirected *)
```

```mathematica
In[1]:= MixedGraphQ[Graph[{1, 2}, {1 -> 2}]]  (* all directed *)
```

### Notes

A graph is *mixed* when it has at least one directed edge and at least one
undirected edge. A purely directed, purely undirected, or edgeless graph is not
mixed.

Several operations — `UndirectedGraph`, `DirectedGraph`, `LineGraph`,
`FindSpanningTree`, the clustering coefficients — are left unevaluated on mixed
graphs, so `MixedGraphQ` is the guard for them.
