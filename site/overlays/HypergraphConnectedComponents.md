### Worked examples

```mathematica
In[1]:= h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
```

```mathematica
In[1]:= Length[HypergraphConnectedComponents[h]]  (* two components: {1..6} and the isolated 7 *)
```

```mathematica
In[1]:= HypergraphConnectedComponents[Hypergraph[{a, b, c, d}, {{a, b}, {c, d}}]]  (* ordered by first vertex *)
```

### Notes

Two vertices are in the same component when a chain of hyperedges joins them;
`HypergraphConnectedComponents` returns these classes as vertex Lists. Isolated
vertices are singleton components.

The ordering is deterministic and differs from Graph `ConnectedComponents`
(largest-first when undirected): hypergraph components are ordered by their first
vertex in `VertexList` order, with vertices inside a component ascending. The
engine is union–find, near-linear in the total incidence. `ConnectedHypergraphQ`
is the companion predicate — one non-empty component.
