### Worked examples

```mathematica
(* ten auto-levels, coloured by height *)
In[1]:= ContourPlot[Sin[x] + Cos[y], {x, -3, 3}, {y, -3, 3}]
```

```mathematica
(* explicit contour values: hyperbolas at the given levels *)
In[1]:= ContourPlot[x^2 - y^2, {x, -2, 2}, {y, -2, 2}, Contours -> {-2, -1, 0, 1, 2}]
```

```mathematica
(* filled contours with a named ramp *)
In[1]:= ContourPlot[x^2 + y^2, {x, -2, 2}, {y, -2, 2}, ContourShading -> True, Contours -> 8]
```

```mathematica
(* an equation plots as its single zero-contour *)
In[1]:= Head[ContourPlot[x^2 + y^2 == 1, {x, -2, 2}, {y, -2, 2}]]
```

```mathematica
(* every level is drawn as marching-squares Line segments *)
In[1]:= Length[Cases[ContourPlot[x^2 + y^2, {x, -2, 2}, {y, -2, 2}, Contours -> 5, PlotPoints -> 25], _Line, Infinity]] > 0
```

### Notes

`ContourPlot` is `HoldAll`. It samples `f` on an `(N+1) x (N+1)` grid (`PlotPoints`
default 75) and runs **marching squares** per contour level, disambiguating the
two saddle cell configurations by the bilinear cell-centre value. Each isoline
segment is emitted as a `Line[...]`. Levels come from `Contours -> n` (a count),
`Contours -> {...}` (explicit values), or the default of 10 interior levels evenly
spaced between the data min and max.

By default (Automatic) only the lines are drawn; shading turns on when a
`ColorFunction` is given, or explicitly with `ContourShading -> True`, filling
each cell with a `Rectangle` coloured by its normalised average.
`ContourPlot[lhs == rhs, ...]` rewrites the body to `lhs - rhs` and draws the
single zero-contour; a list of equations plots each as its own curve. There is no
adaptive refinement, so contour smoothness is set by `PlotPoints`.
