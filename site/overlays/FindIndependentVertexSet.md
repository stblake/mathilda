### Worked examples

```mathematica
In[1]:= FindIndependentVertexSet[PathGraph[{a, b, c, d, e}]]  (* one maximum set, wrapped in a list *)
```

```mathematica
In[1]:= FindIndependentVertexSet[StarGraph[5]]  (* the leaves *)
```

```mathematica
In[1]:= FindIndependentVertexSet[CycleGraph[5]]  (* alpha of a five-cycle is two *)
```

```mathematica
In[1]:= FindIndependentVertexSet[PetersenGraph[]]  (* alpha of the petersen graph is four *)
```

```mathematica
In[1]:= FindIndependentVertexSet[CycleGraph[6], {2}, All]  (* maximal sets of size two, all of them *)
```

### Notes

With only a graph the head returns `{s}`, one maximum independent set. With a size specification (`k`, `{k}` or `{kmin, kmax}`) and an optional count (`n` or `All`) it lists *maximal* independent sets of those sizes, as `FindClique` does.

Which maximum set is returned among ties is deterministic for a given graph. The search is exact but NP-hard in general; past its node/work budget, or a `TimeConstrained` limit, the call stays unevaluated.
