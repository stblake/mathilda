### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

```mathematica
In[1]:= HyperedgeConnectedComponents[h]  (* hyperedge components, as 1-based index lists *)
```

```mathematica
In[1]:= HyperedgeConnectedComponents[Hypergraph[{{1, 2, 3}, {2, 3, 4}, {5, 6}}], 2]  (* with s = 2, joined only on a 2-vertex overlap *)
```

### Notes

Where `HypergraphConnectedComponents` groups vertices, `HyperedgeConnectedComponents`
groups the **hyperedges**: two hyperedges are in the same component when a chain of
pairwise-intersecting hyperedges links them. Components come back as Lists of
1-based hyperedge indices (`EdgeList` positions), ordered by smallest index.

`HyperedgeConnectedComponents[h, s]` gives the s-connected components of Aksoy et
al.: consecutive hyperedges must share at least `s` vertices. Hyperedges with fewer
than `s` distinct vertices lie on no s-walk and are omitted, so empty hyperedges
never appear. `s = 1` is union–find; `s >= 2` counts overlaps. `s` must be a
positive integer.
