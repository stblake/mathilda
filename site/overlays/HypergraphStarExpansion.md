### Worked examples

```mathematica
In[1]:= g = HypergraphStarExpansion[Hypergraph[{{1, 2}, {2, 3}}]]
```

```mathematica
In[1]:= VertexList[g]  (* the original vertices, then one Hyperedge[j] node per hyperedge *)
```

```mathematica
In[1]:= EdgeList[g]  (* each vertex joins the hyperedge nodes it lies in *)
```

### Notes

The star (incidence) expansion is the bipartite `Graph` of the
vertex–hyperedge incidence: its vertices are `VertexList[h]` followed by the nodes
`Hyperedge[1], ..., Hyperedge[m]`, with an edge `v <-> Hyperedge[j]` for each
vertex `v` of hyperedge `j`. Its edge count is the total incidence `Σ|e|`.

The star expansion is also the layout skeleton `HypergraphPlot` draws from. The one
case it declines: if a vertex of `h` happens to *be* a `Hyperedge[j]` node the
construction would be ambiguous, so the call is left unevaluated. It accepts a bare
List of hyperedges as well as a `Hypergraph`.
