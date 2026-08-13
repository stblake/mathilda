# DensityPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DensityPlot[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]`**

Renders f(x,y) as a heatmap: each grid cell is coloured by its function value via ColorFunction (default: thermal blue→yellow ramp). DensityPlot is HoldAll: f is held unevaluated until x and y are bound to numeric values. Returns a Graphics\[...\] object. Options: PlotPoints          grid resolution per axis (default 50) ColorFunction       named ramp string or f\[t\]→color (t in \[0,1\]). Ramps: "Rainbow", "CoolTones", "WarmTones", "Greyscale", "Temperature" (all keyed to normalised z value, t∈\[0,1\]) ColorFunctionScaling True (default): normalise z to \[0,1\] before calling ColorFunction; False: pass raw z RegionFunction      f\[x,y\] mask; excluded cells are not drawn PlotLegends         Automatic: attach a vertical color scale bar Standard Graphics options (Axes, AspectRatio→1, Frame, PlotRange, AxesLabel, GridLines, ImageSize, Background, PlotLabel, …) pass through to the Graphics\[...\] result.

## Examples (4)

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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
