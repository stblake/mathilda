### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

```mathematica
In[1]:= EdgeList[HypergraphLineGraph[h]]  (* an edge joins hyperedges i and j when they intersect *)
```

```mathematica
In[1]:= VertexCount[HypergraphLineGraph[h]]  (* one vertex per hyperedge, even an isolated one *)
```

### Notes

The line graph puts one vertex per hyperedge (`1..m`) and joins `i` and `j` when
hyperedges `i` and `j` intersect. `HypergraphLineGraph[h, s]` gives the **s-line
graph** of Aksoy et al.: `i <-> j` only when the two hyperedges share at least `s`
vertices, the overlap notion that underlies high-order hypergraph walks.

Hyperedges are read as sets, edges are listed in `(i, j)` order, and every
hyperedge index is a vertex even when isolated. The overlaps are counted straight
through the incidence lists — the full pairwise intersection matrix is never
formed — so the cost is the number of hyperedge pairs meeting at a vertex. `s` must
be a positive integer.
