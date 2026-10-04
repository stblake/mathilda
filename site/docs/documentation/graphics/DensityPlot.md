# DensityPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DensityPlot[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]`**

Renders f(x,y) as a heatmap: each grid cell is coloured by its function value via ColorFunction (default: thermal blue→yellow ramp). DensityPlot is HoldAll: f is held unevaluated until x and y are bound to numeric values. Returns a Graphics\[...\] object. Options: PlotPoints          grid resolution per axis (default 50) ColorFunction       named ramp string or f\[t\]→color (t in \[0,1\]). Ramps: "Rainbow", "CoolTones", "WarmTones", "Greyscale", "Temperature" (all keyed to normalised z value, t∈\[0,1\]) ColorFunctionScaling True (default): normalise z to \[0,1\] before calling ColorFunction; False: pass raw z RegionFunction      f\[x,y\] mask; excluded cells are not drawn PlotLegends         Automatic: attach a vertical color scale bar Standard Graphics options (Axes, AspectRatio→1, Frame, PlotRange, AxesLabel, GridLines, ImageSize, Background, PlotLabel, …) pass through to the Graphics\[...\] result.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= DensityPlot[Sin[x] Sin[y], {x, -4, 4}, {y, -3, 3}]
Out[1]= -Graphics-
```

### Options (3)

```mathematica
In[2]:= DensityPlot[x^2 - y^2, {x, -2, 2}, {y, -2, 2}, ColorFunction -> "Rainbow", PlotPoints -> 60]
Out[2]= -Graphics-

In[3]:= DensityPlot[Sin[x + y], {x, 0, 6}, {y, 0, 6}, ColorFunction -> (GrayLevel[#]&), PlotLegends -> Automatic]
Out[3]= -Graphics-

In[4]:= DensityPlot[x^2 + y^2, {x, -3, 3}, {y, -3, 3}, RegionFunction -> Function[{x,y}, x^2 + y^2 < 4]]
Out[4]= -Graphics-
```

### Applications (5)

```mathematica
In[5]:= DensityPlot[Sin[x] Sin[y], {x, -4, 4}, {y, -3, 3}]
Out[5]= -Graphics-

In[6]:= DensityPlot[x^2 - y^2, {x, -2, 2}, {y, -2, 2}, ColorFunction -> "Rainbow", PlotPoints -> 40]
Out[6]= -Graphics-

In[7]:= DensityPlot[x^2 + y^2, {x, -3, 3}, {y, -3, 3}, RegionFunction -> Function[{x, y}, x^2 + y^2 < 4]]
Out[7]= -Graphics-

In[8]:= Head[DensityPlot[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 10]]
Out[8]= Graphics

In[9]:= Attributes[DensityPlot]
Out[9]= {HoldAll, Protected}
```

## Algorithm

densityplot.c — DensityPlot[f, {x,xmin,xmax}, {y,ymin,ymax}, opts...]

Renders f(x,y) as a heatmap: each grid cell is coloured by its average function value via ColorFunction (default: thermal blue→yellow ramp from plot_common's thermal_rgb). Returns Graphics[...] auto-displayed by the REPL.

DensityPlot is HoldAll: f and iterator specs are held unevaluated until x and y get numeric values, matching ContourPlot's semantics.

Options:

```text
  PlotPoints            grid resolution per axis (default 50)
  ColorFunction         f[t] → color, or "Rainbow" / "Temperature"
  ColorFunctionScaling  True (default): normalise z to [0,1]; False: raw z
  RegionFunction        f[x,y] mask; excluded cells are not drawn
  PlotLegends           Automatic: attach a $StreamColorBar colour scale
  Standard Graphics options pass through (Axes, AspectRatio→1, Frame, …) 
```

## Implementation notes

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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/graphics/densityplot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/densityplot.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)

## Notes & additional examples

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
