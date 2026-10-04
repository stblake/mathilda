### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

```mathematica
In[1]:= HypergraphDistance[h, 1, 6]  (* fewest hyperedge hops from vertex 1 to vertex 6 *)
```

```mathematica
In[1]:= HypergraphDistance[h, 1]  (* distances from vertex 1 to every vertex, in VertexList order *)
```

```mathematica
In[1]:= HypergraphDistance[h, 1, 1]  (* the distance from a vertex to itself is 0 *)
```

### Notes

`HypergraphDistance[h, u, v]` is the least number of hyperedges in a chain of
vertices `u = x_0, x_1, ..., x_k = v` in which consecutive vertices share a
hyperedge — equivalently the shortest-path distance in the clique expansion
`HypergraphCliqueExpansion[h]`. It is `0` for `u == v` and `Infinity` when `u` and
`v` lie in different components. `HypergraphDistance[h, u]` gives the whole
distance vector from `u`.

The BFS walks vertices but expands each hyperedge at most once, so it never
materialises the clique graph and stays linear in the total incidence. Contrast
`HyperedgeDistance`, which measures distance *between hyperedges*. An unknown
vertex leaves the call unevaluated.
