### Worked examples

```mathematica
(* overlay two graphics: each input becomes its own directive scope *)
In[1]:= Show[Graphics[{Red, Line[{{0, 0}, {1, 1}}]}], Graphics[Point[{2, 2}]]][[1]]
```

```mathematica
(* combining several plots stays a single Graphics[] object *)
In[1]:= Head[Show[Plot[Sin[x], {x, 0, 1}], Plot[Cos[x], {x, 0, 1}]]]
```

```mathematica
(* PlotRange is the union of the inputs' explicit ranges *)
In[1]:= Cases[Show[Graphics[Line[{{0, 0}, {1, 1}}], PlotRange -> {{0, 1}, {0, 1}}], Graphics[Point[{2, 3}]]], (PlotRange -> v_) :> v]
```

```mathematica
(* merging an option into a single graphics returns the merged object *)
In[1]:= Show[Graphics[{Point[{0, 0}]}], Axes -> True]
```

### Notes

`Show` is the engine's compositor and, via the REPL front end, its display path:
any top-level `Graphics`/`Graphics3D` result is auto-rendered, so `Plot[...]`,
`Show[...]` and `g // Graphics` all reach the window through one route. `Show[g,
opts]` merges options into `g` (later-wins, unknowns appended); `Show[g1, g2,
...]` overlays, wrapping **each input's primitives in its own `List` directive
scope** so a colour or `Dashing` in one input cannot restyle another, and baking
each 2-D input's own `PlotStyle` into its scope.

`PlotRange` of a combination is the **union** of the inputs' ranges — an input
without an explicit range contributes the extent of its primitives — and stays
`Automatic` unless at least one input fixed a range. Options other than the merged
ones are taken from the first graphic. `Show` declines (stays unevaluated) if its
first argument is not a graphics, if a graphics follows an option, or if 2-D and
3-D graphics are mixed.
