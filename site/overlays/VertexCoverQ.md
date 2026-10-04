### Worked examples

```mathematica
In[1]:= VertexCoverQ[CycleGraph[4], {1, 3}]  (* two opposite vertices touch all four edges *)
```

```mathematica
In[1]:= VertexCoverQ[CycleGraph[4], {1, 2}]  (* the opposite edge is left uncovered *)
```

```mathematica
In[1]:= VertexCoverQ[CompleteGraph[3], {1, 2}]  (* every edge of the triangle meets {1,2} *)
```

### Notes

`VertexCoverQ[g, vlist]` is `True` when every edge of `g` has at least one
endpoint in `vlist`. Edge direction is ignored and repeated vertices are
allowed.

It is a total predicate — a non-graph argument or a non-vertex element gives
`False`, never an unevaluated result — and it is the exact complement of
`IndependentVertexSetQ`: `vlist` covers `g` iff the vertices outside `vlist` form
an independent set.
