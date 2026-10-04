### Worked examples

```mathematica
(* a literal grid heatmap on the default Greyscale ramp *)
In[1]:= ArrayPlot[{{1, 0, 1}, {0, 1, 0}, {1, 0, 1}}]
```

```mathematica
(* colour-literal cells are painted directly, so ArrayPlot doubles as a pixel grid *)
In[1]:= ArrayPlot[{{Red, Blue}, {Blue, Red}}]
```

```mathematica
(* ColorRules pin named values, Greyscale elsewhere *)
In[1]:= ArrayPlot[{{1, 0, 0.5}, {0, 1, 0.5}}, ColorRules -> {1 -> Pink, 0 -> Yellow}]
```

```mathematica
(* one Rectangle per cell, with no interpolation between cells *)
In[1]:= Length[Cases[ArrayPlot[{{1, 0}, {0, 1}}], _Rectangle, Infinity]]
```

```mathematica
(* not HoldAll, so the array is evaluated first; Protected *)
In[1]:= Attributes[ArrayPlot]
```

### Notes

`ArrayPlot` renders an already-evaluated array (a nested `List` or an `NDArray`) as
a grid of coloured cells with **no interpolation** — row 1 at the top, column 1 at
the left. It is not `HoldAll`. Each cell is classified as a colour literal
(`RGBColor`/`GrayLevel`/`Hue`/`CMYKColor`, painted as-is) or a plain number (coloured
through `ColorFunction`); the two kinds may mix freely in one array, and a cell that
is neither leaves the call unevaluated.

The default `ColorFunction` is `"Greyscale"` (white at the array minimum, black at
the maximum — matching Mathematica, unlike the Viridis-default sampled plotters).
`ColorRules` overrides exact-matching values before the ramp; `Mesh -> All` draws
grid lines; `PlotLegends -> Automatic` attaches a scale bar, but only when at least
one cell is numeric. `AspectRatio` defaults to `rows/cols`, so cells render square.
