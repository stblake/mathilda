# HypergraphPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphPlot[h, opts] gives a Graphics object drawing the hypergraph h (or a list of hyperedges): vertices placed by the stress layout of the star expansion, each hyperedge of 3 or more vertices a translucent rounded hull in its own palette colour, one of 2 a stadium and one of 1 a circle around its vertex, vertex disks on top. Options: VertexLabels, VertexCoordinates, VertexStyle and GraphLayout as in GraphPlot; others pass through to Graphics.`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Head[HypergraphPlot[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]]]
Out[1]= Graphics

In[2]:= Count[HypergraphPlot[Hypergraph[{{1,2,3},{3,4},{4,5,6},{7}}]], _Polygon, Infinity]
Out[2]= 4

In[3]:= HypergraphPlot[Hypergraph[{{1,2,3}}]] === HypergraphPlot[Hypergraph[{{1,2,3}}]]
Out[3]= True

In[4]:= HypergraphPlot[5]
Out[4]= HypergraphPlot[5]
```

### Options (1)

```mathematica
In[5]:= Count[HypergraphPlot[{{1, 2, 3}, {3, 4}}, VertexLabels -> "Name"], _Text, Infinity]
Out[5]= 4
```

## Implementation notes

- `Protected`. Implemented in `src/graph/hyp_plot.c`; Mathematica has no
  built-in hypergraph drawing (the Function Repository's `HypergraphPlot` is the
  model). Deterministic, like `GraphPlot`.
- **Layout**: the vertices are placed by the stress layout of the star
  expansion (one extra node per hyperedge of two or more vertices, joined to its
  members), so the members of a hyperedge sit around a common centre and
  hyperedges sharing vertices are drawn side by side. `GraphLayout` picks another
  embedding of the star expansion; `VertexCoordinates` overrides any subset.
- **Hyperedges**: each is the convex hull of its (distinct) members, inflated
  by a margin with rounded corners (sampled every 15 degrees), drawn as a
  translucent (`Opacity[0.22]`) filled `Polygon` with a darker outline, in its
  own colour of the `ColorData[97]` palette (cycled by hyperedge index). A
  hyperedge of size 2 is therefore a stadium (a thick translucent line with round
  caps) and one of size 1 a circle around its vertex; an empty one is not drawn.
  Larger shapes are drawn first so smaller ones stay visible; a hyperedge sharing
  a vertex with smaller ones gets a wider margin, so nested and repeated
  hyperedges show as concentric outlines.
- **Vertices** are dark disks drawn on top; labels avoid the directions of the
  vertex's hyperedges.
- Unevaluated on a non-hypergraph.

**Attributes:** `Protected`.

## References

**See also:** [GraphPlot](../../graphs/GraphPlot/)

- Source: [`src/graph/hyp_init.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_init.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_graphplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphplot.c)
