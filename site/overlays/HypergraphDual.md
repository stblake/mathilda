### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

```mathematica
In[1]:= InputForm[HypergraphDual[h]]  (* vertices 1..m, one hyperedge per vertex of h *)
```

```mathematica
In[1]:= {VertexCount[HypergraphDual[h]], EdgeCount[HypergraphDual[h]]}  (* the dual has m vertices and n hyperedges *)
```

### Notes

The dual exchanges the roles of vertices and hyperedges. Its vertices are
`1, ..., m` (one per hyperedge of `h`), and for each vertex `v` of `h`, in
`VertexList` order, it has the hyperedge listing the ascending indices of the
hyperedges of `h` that contain `v`. An isolated vertex of `h` becomes an empty
hyperedge of the dual. Hence the dual of a hypergraph with `n` vertices and `m`
hyperedges has `m` vertices and `n` hyperedges.

Hyperedges are read as sets, so `HypergraphDual[HypergraphDual[h]]` recovers the
incidence structure of `h` on `1..n`. The build is linear in the total incidence.
A bare List of hyperedges is not accepted — the head needs the memoized incidence
of a real `Hypergraph` object.
