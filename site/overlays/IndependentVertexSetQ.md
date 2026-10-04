### Worked examples

```mathematica
In[1]:= IndependentVertexSetQ[CycleGraph[4], {1, 3}]  (* opposite corners of a 4-cycle are non-adjacent *)
```

```mathematica
In[1]:= IndependentVertexSetQ[CycleGraph[4], {1, 2}]  (* adjacent vertices -- not independent *)
```

```mathematica
In[1]:= IndependentVertexSetQ[CompleteGraph[4], {2}]  (* a single vertex is always independent *)
```

### Notes

`IndependentVertexSetQ[g, vlist]` is `True` when no two vertices of `vlist` are
adjacent in `g`. Edge direction is ignored, so an `UndirectedEdge` and a
`DirectedEdge` both count as an adjacency, and repeated vertices in `vlist` are
allowed.

It is a total predicate: a non-graph first argument, or a `vlist` element that is
not a vertex of `g`, gives `False` rather than staying unevaluated. It is the
complement of `VertexCoverQ` — `vlist` is independent iff its complement is a
vertex cover.
