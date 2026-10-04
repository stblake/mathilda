### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

```mathematica
In[1]:= InputForm[HypergraphVertexDelete[h, 3]]  (* also removes every hyperedge containing 3 *)
```

```mathematica
In[1]:= InputForm[HypergraphVertexDelete[h, {1, 7}]]  (* delete several vertices at once *)
```

### Notes

`HypergraphVertexDelete[h, v]` removes the vertex `v` **and every hyperedge
containing it**, exactly as `VertexDelete` does for a `Graph`; the list form
removes several. Deleting vertex `3` from the running example therefore drops both
hyperedges `{1, 2, 3}` and `{3, 4}`.

A named vertex that is not in `h` leaves the call unevaluated. As with
`HypergraphVertexAdd`, any `List` argument is a list of vertices, and the edit
produces a fresh object without mutating the original. To delete hyperedges while
keeping their vertices, use `HypergraphEdgeDelete`; to intersect hyperedges with a
vertex set instead of dropping them, use `HypergraphRestriction`.
