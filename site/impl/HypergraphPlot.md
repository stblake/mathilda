---
references:
  - "A. M. Andrew, *Another efficient algorithm for convex hulls in two dimensions*, Information Processing Letters **9**(5) (1979) 216-219 (monotone chain)."
source: src/graph/hyp_plot.c
---
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
