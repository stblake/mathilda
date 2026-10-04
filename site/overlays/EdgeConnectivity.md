### Worked examples

```mathematica
In[1]:= EdgeConnectivity[CycleGraph[5]]  (* two edges must go to break a cycle *)
```

```mathematica
In[1]:= EdgeConnectivity[CompleteGraph[4]]  (* equals the minimum degree, 3 *)
```

```mathematica
In[1]:= EdgeConnectivity[PathGraph[{1, 2, 3}]]  (* a single bridge disconnects a path *)
```

```mathematica
In[1]:= EdgeConnectivity[CompleteGraph[4], 1, 2]  (* the minimum 1-2 edge cut *)
```

### Notes

`EdgeConnectivity[g]` is the weight of a global minimum edge cut — the cheapest
set of edges whose removal disconnects the graph (strongly, for a directed
graph). `EdgeConnectivity[g, s, t]` restricts this to cuts that separate `s`
from `t`, which by the max-flow/min-cut theorem equals the maximum `s-t` flow
under unit (or `EdgeWeight`) capacities.

With integer weights the answer is an exact integer. A graph with fewer than
two vertices is left unevaluated.
