---
source: src/graphics/densityplot.c
---
**Algorithm.** `builtin_densityplot` is `HoldAll`. It samples `f` on an
`(N+1) x (N+1)` grid (`N = PlotPoints`, default **50**), compiling the body once
with `autocompile_new(body, {x, y}, 2)` and falling back to the interpreter per
sample, then shades each cell by the **average of its four corners**
`(v00+v10+v11+v01)/4` and emits it as a `Rectangle`. Sampling is uniform in
scaled (`ScalingFunctions`) world space. The colour comes from `dp_color`: with
no `ColorFunction` the default is the `default_ramp_rgb` (Viridis) ramp keyed to
the normalised value; a named string resolves through `named_color_ramp`; a user
function is called at arity 1 (`f[t]`). `ColorFunctionScaling -> True` (default)
maps the value to `[0,1]` as `(avg - zmin)/zspan` clamped; `False` passes the raw
value unclamped (clamping a raw value would silently flatten out-of-`[0,1]` data).
`RegionFunction` is tested at the cell centre in data coordinates and excluded
cells are not drawn; each `Rectangle` overlaps `2*du` into its later-drawn
neighbours to close sub-pixel seams. `PlotLegends -> Automatic` appends a
`$StreamColorBar[zmin, zmax, cfn]` metadata node for the renderer's colour scale.
The result is an inert `Graphics[prims, opts...]`.

**Data structures.** `DensityOpts`; a flat `double* grid` of `(N+1)^2` samples;
`Expr** prims` sized `N^2*2 + 4` (a colour directive and a `Rectangle` per cell).
Accepted colour heads are `RGBColor`/`GrayLevel`/`Hue`/`CMYKColor`.

**Complexity / limits.** `O((N+1)^2)` evaluations (~2.6k samples at the default
50); no adaptive sampling, so resolution is set by `PlotPoints`. Defaults
`Frame -> True`, `Axes -> False`, `AspectRatio -> 1`. Unlike `ArrayPlot`, the grid
is a sampled continuous function, not a literal array.
