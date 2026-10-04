# HypergraphPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`HypergraphPlot[h, opts] gives a Graphics object drawing the hypergraph h (or a list of hyperedges): vertices placed by the stress layout of the star expansion, each hyperedge of 3 or more vertices a translucent rounded hull in its own palette colour, one of 2 a stadium and one of 1 a circle around its vertex, vertex disks on top. Options: VertexLabels, VertexCoordinates, VertexStyle and GraphLayout as in GraphPlot; others pass through to Graphics.`**

## Examples (8)

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

### Applications (3)

The result is a Graphics object

```mathematica
In[6]:= Head[HypergraphPlot[Hypergraph[{{1, 2, 3}, {3, 4}}]]]
Out[6]= Graphics
```

One filled hull per hyperedge

```mathematica
In[7]:= Count[HypergraphPlot[Hypergraph[{{1, 2, 3}, {3, 4}, {7}}]], _Polygon, Infinity]
Out[7]= 3
```

Deterministic, like GraphPlot

```mathematica
In[8]:= HypergraphPlot[Hypergraph[{{1, 2}}]] === HypergraphPlot[Hypergraph[{{1, 2}}]]
Out[8]= True
```

## Implementation notes

**Algorithm.** `builtin_hypergraph_plot` returns a `Graphics[...]` object,
deterministically (like `GraphPlot`). Vertices are placed by the **stress layout**
(`glayout.c`) of the star expansion — one extra node per hyperedge of two or more
vertices joined to its members — so a hyperedge's members cluster around a common
centre; `GraphLayout` selects another embedding and `VertexCoordinates` overrides
any subset. Each hyperedge is drawn as the **convex hull** of its members (Andrew's
monotone chain, `convex_hull`) inflated by a margin with rounded corners (the
Minkowski sum with a disk, sampled every 15°, `rounded_hull`): a translucent
(`Opacity 0.22`) filled `Polygon` plus a darker outline in one colour of
`ColorData[97]`, cycled by hyperedge index. Size 2 becomes a stadium, size 1 a
circle. Shapes are drawn largest-first (by inflated area) and a hyperedge sharing a
vertex with smaller ones gets a wider margin (`level`), so nested hyperedges read
as concentric outlines. Vertex disks and labels (direction-avoiding) go on top.

**Data structures.** The `HypView` distinct sets; the star-expansion edge arrays
`eu/ev` and the layout coordinate buffer `xy`; per-hyperedge `pts`/`hull` point
buffers and an `HOrder` draw-order table; the shared graphics-primitive builder
`GDPrims`. The option plumbing (`VertexLabels`, `VertexStyle`, `VertexCoordinates`,
`GraphLayout`) reuses `GraphPlot`'s helpers.

**Complexity / limits.** Dominated by the stress layout; hull and inflation are
`O(k log k)` per hyperedge. A non-hypergraph (and a non-List, non-rule option)
leaves the call unevaluated. Mathematica has no built-in hypergraph drawing; the
Function Repository's `HypergraphPlot` is the model.

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

- A. M. Andrew, *Another efficient algorithm for convex hulls in two dimensions*, Information Processing Letters **9**(5) (1979) 216-219 (monotone chain).
- Source: [`src/graph/hyp_plot.c`](https://github.com/stblake/mathilda/blob/main/src/graph/hyp_plot.c)
- Specification: [`docs/spec/builtins/hypergraphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/hypergraphs.md)
- Tests: [`tests/test_graphplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphplot.c)

## Notes & additional examples

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
