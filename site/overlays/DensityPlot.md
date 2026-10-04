### Worked examples

```mathematica
(* a heatmap of f(x, y), one shaded cell per grid square *)
In[1]:= DensityPlot[Sin[x] Sin[y], {x, -4, 4}, {y, -3, 3}]
```

```mathematica
(* a named colour ramp at higher resolution *)
In[1]:= DensityPlot[x^2 - y^2, {x, -2, 2}, {y, -2, 2}, ColorFunction -> "Rainbow", PlotPoints -> 40]
```

```mathematica
(* a RegionFunction mask leaves excluded cells undrawn *)
In[1]:= DensityPlot[x^2 + y^2, {x, -3, 3}, {y, -3, 3}, RegionFunction -> Function[{x, y}, x^2 + y^2 < 4]]
```

```mathematica
(* the result is an inert Graphics[] object *)
In[1]:= Head[DensityPlot[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 10]]
```

```mathematica
(* HoldAll and Protected *)
In[1]:= Attributes[DensityPlot]
```

### Notes

`DensityPlot` is `HoldAll`. It samples `f` on an `(N+1) x (N+1)` grid (`PlotPoints`
default 50) and shades each cell as a `Rectangle` coloured by the **average of its
four corners** — a sampled continuous function, as opposed to `ArrayPlot`'s literal
array with no interpolation. The default `ColorFunction` is the perceptually-uniform
Viridis ramp.

`ColorFunctionScaling -> True` (the default) normalises each cell value to `[0,1]`
before colouring; `ColorFunctionScaling -> False` passes the raw value unclamped, so
the colour function owns its own domain — which is how a shared colour scale across
several frames is built. `PlotLegends -> Automatic` attaches a vertical colour scale
bar.
