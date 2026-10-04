### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

```mathematica
In[1]:= EdgeList[HypergraphCliqueExpansion[h]]  (* an edge for every pair of vertices sharing a hyperedge *)
```

```mathematica
In[1]:= VertexCount[HypergraphCliqueExpansion[h]]  (* the 2-section keeps every vertex *)
```

### Notes

The clique expansion (2-section) replaces each hyperedge by a clique: `u <-> v` in
the resulting `Graph` exactly when `u` and `v` lie together in some hyperedge. A
vertex repeated inside a hyperedge produces no self-loop, since the hyperedge is
read as its distinct-vertex set, and every vertex of `h` is kept even if isolated.

This is precisely the graph whose shortest-path distance `HypergraphDistance`
measures. Edges are listed in order of first co-occurrence. The head also accepts a
bare List of hyperedges, as the Function Repository function does.
