### Worked examples

```mathematica
In[1]:= PageRankCentrality[CycleGraph[4]]  (* a symmetric graph gives equal ranks *)
```

```mathematica
In[1]:= PageRankCentrality[StarGraph[4]]  (* the hub collects the most rank at the default damping 0.85 *)
```

```mathematica
In[1]:= PageRankCentrality[StarGraph[4], 1/2]  (* lower damping flattens the ranks toward uniform *)
```

```mathematica
In[1]:= PageRankCentrality[Graph[{1, 2, 3}, {DirectedEdge[1, 2], DirectedEdge[2, 3]}]]  (* the sink at the end of a directed path ranks highest *)
```

```mathematica
In[1]:= Total[PageRankCentrality[PathGraph[{1, 2, 3, 4, 5}]]]  (* ranks always sum to one *)
```

### Notes

The optional second argument is the damping factor `a` with `0 <= a <= 1`; the default is `0.85`. The ranks solve `x = a P^T x + (1 - a)/n` and sum to 1. A vertex with no outgoing edge redistributes its rank uniformly over all vertices.

Iteration runs to an L1 change below `10^-14`, tighter than Mathematica's stopping rule, so the two agree to about nine digits. Edge weights are ignored.
