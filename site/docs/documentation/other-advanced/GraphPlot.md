# GraphPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphPlot[g] gives a Graphics object drawing the graph g with a circular vertex layout. Vertex labels are off by default; pass VertexLabels -> True to draw them (in black).`**

## Examples

_No verified examples yet for this function._

## Algorithm

graphplot.c - GraphPlot[g]: render a graph as a Graphics[...] expression.

Emits the same primitives the plotting engine uses (Line, Disk, Text), so the existing renderer draws it with no renderer changes (and the text placeholder is used when USE_GRAPHICS=0). Vertices are laid out on a circle (MVP layout; a force-directed spring layout is the documented future hook). Each edge is a Line between its endpoints, each vertex a Disk plus a Text label.

Directed edges are drawn as plain lines in the MVP (no arrowheads yet).

Memory (SPEC section 4): returns a freshly-allocated Graphics tree; the evaluator frees res.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
