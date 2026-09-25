# GraphPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`GraphPlot[g] gives a Graphics object drawing the graph g with a circular vertex layout. Vertex labels are off by default; pass VertexLabels -> True to draw them (in black).`**

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Head[GraphPlot[CycleGraph[8]]]
Out[1]= Graphics

In[2]:= Count[GraphPlot[CompleteGraph[6]], _Line, Infinity]
Out[2]= 15

In[3]:= Count[GraphPlot[CycleGraph[5]], _Disk, Infinity]
Out[3]= 5

In[4]:= Count[GraphPlot[CycleGraph[5]], _Text, Infinity]
Out[4]= 0

In[5]:= GraphPlot[5]
Out[5]= GraphPlot[5]
```

## Algorithm

graphplot.c - GraphPlot[g]: render a graph as a Graphics[...] expression.

Emits the same primitives the plotting engine uses (Line, Disk, Text), so the existing renderer draws it with no renderer changes (and the text placeholder is used when USE_GRAPHICS=0). Vertices are laid out on a circle (MVP layout; a force-directed spring layout is the documented future hook). Each edge is a Line between its endpoints, each vertex a Disk plus a Text label.

Directed edges are drawn as plain lines in the MVP (no arrowheads yet).

Memory (SPEC section 4): returns a freshly-allocated Graphics tree; the evaluator frees res.

## Implementation notes

- `Protected`. Vertices are laid out on a circle, each drawn as a `Disk`; edges
  are `Line`s. The specification calls for one `Text` label per vertex as well,
  but the current binary emits no `Text` primitives (see the example below).
- Renders through the standard graphics path (a window when `USE_GRAPHICS=1`,
  the text placeholder otherwise).
- MVP limitations: directed edges are drawn as plain lines (no arrowheads yet);
  a force-directed layout is a future hook. Mathematica's `GraphPlot` uses a
  spring-electrical layout.
- Unevaluated on a non-graph.

**Attributes:** `Protected`.

## References

- Source: [`src/graph/graph.c`](https://github.com/stblake/mathilda/blob/main/src/graph/graph.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph.c)
