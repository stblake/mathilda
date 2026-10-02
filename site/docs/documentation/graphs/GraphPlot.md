# GraphPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphPlot[g, opts] gives a Graphics object drawing the graph g (or a list of rules {u -> v, ...}). The default layout is deterministic: tidy layered trees for branching forests, layered drawings for DAGs, stress majorization otherwise, with components packed side by side. Options: GraphLayout -> "StressEmbedding" | "SpringElectricalEmbedding" | "CircularEmbedding" | "LayeredEmbedding" | "BipartiteEmbedding" | "GridEmbedding"; VertexCoordinates -> {{x,y}, ...} or {v -> {x,y}, ...}; VertexLabels -> None | "Name" | Automatic | {v -> lbl, ...}; GraphHighlight -> {v, e, ...} (red, thicker); VertexStyle and EdgeStyle -> a colour or {item -> colour, ...}; EdgeLabels -> "EdgeWeight" | {e -> lbl}; VertexSize -> d (diameter in edge lengths). Directed edges get arrowheads that stop at the target vertex. Other options (ImageSize, PlotLabel, ...) pass through to Graphics.`**

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= Head[GraphPlot[PetersenGraph[]]]
Out[1]= Graphics

In[2]:= Count[GraphPlot[CompleteGraph[6]], _Line, Infinity]
Out[2]= 15

In[3]:= Count[GraphPlot[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]], _Arrow, Infinity]
Out[3]= 3

In[4]:= {Axes, AspectRatio} /. Rest[List @@ GraphPlot[CycleGraph[4]]]
Out[4]= {False, Automatic}

In[5]:= Length[Union[Cases[GraphPlot[GridGraph[{4, 4}]], Disk[{x_, _}, _] :> Round[x, 0.001], Infinity]]]
Out[5]= 4

In[6]:= GraphPlot[5]
Out[6]= GraphPlot[5]
```

### Options (3)

```mathematica
In[7]:= Cases[GraphPlot[PathGraph[{1, 2, 3}], VertexCoordinates -> {{0, 0}, {1, 0}, {2, 1}}], Disk[p_, _] :> p, Infinity]
Out[7]= {{0.0, 0.0}, {1.0, 0.0}, {2.0, 1.0}}

In[8]:= Count[GraphPlot[CycleGraph[5], VertexLabels -> "Name"], _Text, Infinity]
Out[8]= 5

In[9]:= MemberQ[GraphPlot[CycleGraph[3], GraphHighlight -> {1}], RGBColor[1., 0., 0.], Infinity]
Out[9]= True
```

## Algorithm

graphplot.c - GraphPlot[g, opts]: draw a graph as a Graphics[...] expression, plus the drawing helpers GraphPlot shares with HypergraphPlot (glayout.h).

### Pipeline

```text
  1. Layout (glayout.c): GraphLayout -> Automatic (tidy trees for branching
     forests, layered DAGs, stress majorization otherwise), or an explicit
     embedding; VertexCoordinates overrides any subset of the positions.
  2. Size: vertex disks get a radius scaled to the median edge length and the
     layout extent, so a 10-vertex and a 1000-vertex graph both read.
  3. Frame (gd_frame): labels are placed beside their vertex on the side
     with the widest angular gap between incident edges, and the world box
     and page size are solved together so labels are never clipped.
  4. Emit: edges (Line, or Arrow shortened to stop at the target disk, with
     Arrowheads sized to the disks; mutual pairs u->v, v->u offset apart),
     then vertex disks with a thin darker rim, then labels.
```

The output uses only primitives every renderer draws -- Line, Arrow, Disk, Circle, Text, colour and Thickness/Arrowheads directives -- with AspectRatio -> Automatic, Axes -> False and an explicit PlotRange/ImageSize. Colours: vertices RGBColor[0.368417, 0.506779, 0.709798] (ColorData[97][1]), edges a medium grey-blue, GraphHighlight in red and thicker, as Mathematica.

Memory (SPEC section 4): returns a freshly-allocated Graphics tree; the evaluator frees res. Nothing borrowed from res outlives the call.

## Implementation notes

- `Protected`. Implemented in `src/graph/graphplot.c` over the layout engine
  `src/graph/glayout.c`. **Deterministic**: no random numbers anywhere, so the
  same graph always gives an identical `Graphics` expression.
- **Default layout** (`GraphLayout -> Automatic`): a forest with a branching
  vertex is drawn as tidy layered trees (hanging from the tree centre, or from
  the source of an arborescence); an all-directed acyclic graph as a layered
  drawing; everything else (paths and cycles included) by **stress
  majorization** (SMACOF on BFS graph distances, weights `d^-2`), started from
  Pivot MDS. Components of up to 60 vertices also try circle starts and Tutte
  (barycentric) starts from shortest cycles, and keep, among the drawings within
  10% of the least stress, the one with the fewest edge crossings; a
  crossing-free drawing may cost up to 2.5x the stress when it removes at least
  four crossings (so the dodecahedron comes out as its Schlegel diagram while
  the cube stays the textbook Necker cube). Each drawing is rotated to a
  canonical orientation (principal axis horizontal; snapped to the axes when the
  edges are four-fold, then straightened so grids come out exactly as grids;
  vertex 1 on top when there is no preferred axis). Every connected component is
  laid out on its own and the components are shelf-packed, largest first, with
  isolated vertices gathered into a square block.
- **Cost**: full stress majorization up to 1000 vertices per component, Pivot
  MDS (50 pivots, `O(k (n + m))`) above that. A 500-vertex random graph takes
  about 0.15 s, a 25x20 grid 0.03 s.
- `GraphLayout` values: `"StressEmbedding"`, `"SpringElectricalEmbedding"` (Hu's
  spring-electrical model, exact repulsion up to 1000 vertices, grid
  cut-off above), `"CircularEmbedding"` (`VertexList` order, vertex 1 on top),
  `"LayeredEmbedding"` / `"LayeredDigraphEmbedding"` (longest-path layers for a
  DAG, BFS layers from the centre otherwise, dummy vertices on long edges,
  barycentre crossing reduction, isotonic-regression x placement; wide shallow
  drawings get taller layer spacing), `"BipartiteEmbedding"` (the two parts in
  two columns, barycentre-ordered; falls back to stress for a non-bipartite
  graph), `"GridEmbedding"` (`VertexList` order on a square grid).
- `VertexCoordinates -> {{x1, y1}, ...}` (one pair per vertex, `VertexList`
  order) fixes the drawing; `VertexCoordinates -> {v -> {x, y}, ...}` fixes the
  given vertices and lays out the rest.
- `VertexLabels -> None` (default) | `"Name"` | `Automatic` | `True` | `All`
  labels each vertex with its name; `{v -> lbl, ...}` labels only those. A label
  sits beside its vertex on the side with the widest angular gap between the
  incident edges (upper right when free), in 10 pt Helvetica, and the frame
  grows so no label is clipped.
- **Directed edges** are `Arrow`s with an `Arrowheads` directive sized to the
  vertex disks; each arrow starts outside its source disk and its tip stops just
  short of the target disk. A mutual pair `u -> v`, `v -> u` is drawn as two
  arrows offset to either side.
- `GraphHighlight -> {v, ..., e, ...}`: highlighted vertices and edges are drawn
  red (`RGBColor[1, 0, 0]`), vertices 15% larger and edges 2.5x thicker, on top
  of the others. Edges may be written `u <-> v`, `UndirectedEdge[u, v]`,
  `u -> v` or `DirectedEdge[u, v]`.
- `VertexStyle -> style` or `{v -> style, ...}`; `EdgeStyle -> style` or
  `{e -> style, ...}` (a colour, or a list of directives), e.g. a vertex
  colouring from `FindVertexColoring`.
- `EdgeLabels -> "EdgeWeight"` writes each weight at its edge midpoint (offset
  towards the inside of the drawing); `EdgeLabels -> {e -> lbl, ...}` labels
  chosen edges. `VertexSize -> d` sets the disk diameter to `d` edge lengths.
- **Look**: vertices are disks in `RGBColor[0.368417, 0.506779, 0.709798]`
  (Mathematica's `ColorData[97]` blue) with a thin darker rim, edges 1.1 pt in
  grey-blue `RGBColor[0.571589, 0.586483, 0.699215]`, the disk radius scaled to
  the median edge length and the extent of the drawing. The result carries
  `PlotRange`, `AspectRatio -> Automatic`, `Axes -> False` and an explicit
  `ImageSize -> {w, h}` (300 pt on the long side, growing gently with the vertex
  count), so `Export["g.pdf", GraphPlot[g]]` gives a tight, equal-aspect picture.
- Unevaluated on a non-graph, or when an argument after the graph is not a rule.

**Attributes:** `Protected`.

## References

**See also:** [ImageSize](../../other-advanced/ImageSize/), [VertexList](../../graphs/VertexList/), [FindVertexColoring](../../graphs/FindVertexColoring/)

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
- Tests: [`tests/test_graphplot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphplot.c)
