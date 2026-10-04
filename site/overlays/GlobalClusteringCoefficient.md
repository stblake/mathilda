### Worked examples

```mathematica
In[1]:= GlobalClusteringCoefficient[CompleteGraph[4]]  (* every triple closes, so 1 *)
```

```mathematica
In[1]:= GlobalClusteringCoefficient[CycleGraph[5]]  (* no triangles, so 0 *)
```

```mathematica
In[1]:= GlobalClusteringCoefficient[PetersenGraph[]]  (* the Petersen graph is triangle-free *)
```

### Notes

The global (or *transitivity*) coefficient is `3 x (triangles) /
(connected triples)`: the fraction of paths of length two that are closed into
a triangle. It is `1` for a complete graph and `0` for any triangle-free graph.

It is computed exactly as a `Rational`. This is a single graph-wide ratio and
is in general not equal to `MeanClusteringCoefficient`, which averages the
per-vertex local coefficients.
