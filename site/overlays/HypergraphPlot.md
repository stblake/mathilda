### Worked examples

```mathematica
In[1]:= Head[HypergraphPlot[Hypergraph[{{1, 2, 3}, {3, 4}}]]]  (* the result is a Graphics object *)
```

```mathematica
In[1]:= Count[HypergraphPlot[Hypergraph[{{1, 2, 3}, {3, 4}, {7}}]], _Polygon, Infinity]  (* one filled hull per hyperedge *)
```

```mathematica
In[1]:= HypergraphPlot[Hypergraph[{{1, 2}}]] === HypergraphPlot[Hypergraph[{{1, 2}}]]  (* deterministic, like GraphPlot *)
```

### Notes

`HypergraphPlot[h]` returns a `Graphics[...]` object. Vertices are placed by the
stress layout of the star expansion, so a hyperedge's members cluster around a
shared centre and hyperedges that share vertices are drawn side by side. Each
hyperedge is the convex hull of its members, inflated with rounded corners and
drawn as a translucent filled `Polygon` with a darker outline, one palette colour
per hyperedge; a hyperedge of two vertices becomes a stadium and one of a single
vertex a circle. Larger shapes are drawn first and overlapping hyperedges get
widening margins, so nested ones read as concentric outlines.

The examples above query the result (`Head`, `Count`, equality) rather than
displaying it, so they are safe to run headless; evaluating `HypergraphPlot[h]`
directly renders the picture. Options `VertexLabels`, `VertexCoordinates`,
`VertexStyle` and `GraphLayout` behave as in `GraphPlot`; others pass through to
`Graphics`. Output is deterministic. `HypergraphPlot[{e1, e2, ...}]` draws a bare
list of hyperedges.
