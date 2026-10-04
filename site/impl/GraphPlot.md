---
references:
  - "E. R. Gansner, Y. Koren and S. North, *Graph drawing by stress majorization*, in Graph Drawing (GD 2004), LNCS **3383**, Springer (2005) 239-250."
  - "Y. Hu, *Efficient, high-quality force-directed graph drawing*, The Mathematica Journal **10** (2005) 37-71."
  - "T. M. J. Fruchterman and E. M. Reingold, *Graph drawing by force-directed placement*, Software: Practice and Experience **21** (1991) 1129-1164."
source: src/graph/graphplot.c
---
**Algorithm.** `builtin_graph_plot` draws a graph as a `Graphics[...]` expression
through a four-stage pipeline. **(1) Layout** (`glayout_compute`, in
`src/graph/glayout.c`): `GraphLayout -> Automatic` chooses by structure — a tidy
tree for a branching forest (leaves in DFS order, parents centred over their
children), a layered drawing for a DAG, and **stress majorization** (SMACOF on
BFS graph distances) for everything else; the explicit methods are
`"CircularEmbedding"`, `"SpringElectricalEmbedding"` (Hu's spring-electrical model
with Fruchterman-Reingold grid repulsion), `"StressEmbedding"`,
`"LayeredEmbedding"`/`"LayeredDigraphEmbedding"`, `"BipartiteEmbedding"` and
`"GridEmbedding"`, an inapplicable one falling back to stress. `VertexCoordinates`
overrides any subset of the positions, and a complete spec skips the layout
entirely. **(2) Size**: the vertex-disk radius is scaled to the median edge length
and the layout extent, so a 10- and a 1000-vertex graph both read. **(3) Frame**
(`gd_frame`): each label is placed beside its vertex on the side with the widest
angular gap between incident edges, and the world box and page size are solved
together so labels are never clipped. **(4) Emit**: edges first (a `Line` when
undirected, an `Arrow` shortened to stop at the target disk with `Arrowheads`
sized to the disks, mutual pairs `u->v`/`v->u` offset apart), then vertex `Disk`s
with a thin darker rim, then labels. `GraphPlot[{rules}]` plots `Graph[{rules}]`;
options the head does not consume pass through to `Graphics`.

**Data structures.** Edges are read via `graph_edge_indices` into integer endpoint
arrays (`eu`/`ev`) and a per-edge direction-bit array; coordinates live in a flat
`double xy[]` buffer filled by the layout. Primitives accumulate in a `GDPrims`
vector and are finished (`gd_finish`) into a `Graphics[...]` tree built from only
the primitives every renderer draws — `Line`, `Arrow`, `Disk`, `Circle`, `Text`
and colour/`Thickness`/`Arrowheads` directives — with `AspectRatio -> Automatic`,
`Axes -> False` and an explicit `PlotRange`/`ImageSize`. Vertex colour is
`RGBColor[0.368417, 0.506779, 0.709798]` (`ColorData[97][1]`), edges a medium
grey-blue, `GraphHighlight` red and thicker, matching Mathematica.

**Complexity / limits.** Cost is dominated by the layout: stress majorization is
`O(n^2)` per sweep (large components skip the per-sweep majorization), the
all-pairs BFS distance matrix is `O(n (V + E))`, and the exact graph centre used
to seed BFS layouts is computed only up to 2000 vertices. The output is a
`Graphics` object — a renderable drawing, not a machine value — so, like other
graphics heads, `GraphPlot` carries no NDArray/`Compile` fast path.
