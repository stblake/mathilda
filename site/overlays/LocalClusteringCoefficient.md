### Worked examples

```mathematica
In[1]:= LocalClusteringCoefficient[CompleteGraph[4]]  (* every neighbourhood is complete, so all 1 *)
```

```mathematica
In[1]:= LocalClusteringCoefficient[CycleGraph[5]]  (* no neighbour pair is adjacent *)
```

```mathematica
In[1]:= LocalClusteringCoefficient[WheelGraph[5], 1]  (* the hub: 4 of 6 neighbour pairs are joined *)
```

### Notes

The local clustering coefficient of a vertex is how close its neighbourhood is
to a clique: the number of edges among its neighbours divided by the number of
possible such edges, `C(d, 2)`. A vertex of degree below 2 gets `0`.

`LocalClusteringCoefficient[g]` returns one exact `Rational` per vertex in
`VertexList` order; `LocalClusteringCoefficient[g, v]` gives the single value.
The mean of these numbers is `MeanClusteringCoefficient`.
