### Worked examples

```mathematica
In[1]:= Head[GraphPlot[CycleGraph[4]]]  (* the result is an ordinary Graphics object *)
```

```mathematica
In[1]:= Head[GraphPlot[{1 -> 2, 2 -> 3, 3 -> 1}]]  (* a list of rules is plotted as its Graph *)
```

```mathematica
In[1]:= MatchQ[GraphPlot[CompleteGraph[4], GraphLayout -> "CircularEmbedding"], _Graphics]  (* an explicit layout method *)
```

```mathematica
In[1]:= Head[GraphPlot[PathGraph[{1, 2, 3}], VertexLabels -> "Name"]]  (* labels sit beside each vertex on its widest gap *)
```

### Notes

`GraphPlot[g]` returns a `Graphics[...]` expression drawing `g`, built from plain
primitives (`Line`, `Arrow`, `Disk`, `Circle`, `Text`), so it renders anywhere
`Graphics` does. Because that object is large and layout-dependent, the examples
above probe it with `Head` or `MatchQ` rather than printing the drawing itself.

The vertex positions come from `GraphLayout`: `Automatic` picks a tidy tree for a
forest, a layered drawing for a DAG, and stress-majorization placement otherwise,
while `"CircularEmbedding"`, `"SpringElectricalEmbedding"`, `"StressEmbedding"`,
`"LayeredEmbedding"`, `"BipartiteEmbedding"` and `"GridEmbedding"` request a
specific method. `VertexCoordinates` pins any subset of the vertices, and options
`GraphPlot` does not itself consume (such as `ImageSize` or `PlotLabel`) pass
through to the enclosing `Graphics`. A bare list of edge rules is plotted as the
corresponding `Graph`.
