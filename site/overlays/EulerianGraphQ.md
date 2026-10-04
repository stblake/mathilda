### Worked examples

```mathematica
In[1]:= EulerianGraphQ[CycleGraph[5]]  (* every degree is 2, so a cycle is Eulerian *)
```

```mathematica
In[1]:= EulerianGraphQ[PathGraph[{1, 2, 3}]]  (* the two ends have odd degree *)
```

```mathematica
In[1]:= EulerianGraphQ[CompleteGraph[5]]  (* every vertex has even degree 4 *)
```

```mathematica
In[1]:= EulerianGraphQ[CompleteGraph[4]]  (* degree 3 is odd, so not Eulerian *)
```

### Notes

A graph is Eulerian when it has a closed walk traversing every edge once. The
condition is: all edges in one component, and every vertex of even degree
(undirected) or with equal in- and out-degree (directed). `K_n` is Eulerian
exactly when `n` is odd.

Mixed graphs (both directed and undirected edges) are left unevaluated. Use
`FindEulerianCycle` to obtain an actual tour.
