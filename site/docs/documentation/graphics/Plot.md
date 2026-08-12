# Plot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Plot[f, {x, xmin, xmax}, opts...]`**

Adaptively samples f over \[xmin, xmax\], displays the resulting curve in an interactive window, and returns it as a Graphics\[...\] object. A list of functions Plot\[{f1, f2, ...}, {x, xmin, xmax}\] draws each on the same axes in a distinct palette colour. Options: PlotPoints (initial sample count, default 50), MaxRecursion (adaptive refinement depth, default 6), MaxPlotPoints (overall point cap, default Infinity), Mesh (All overlays the evaluation points as dots; default None), PlotRange, PlotRangePadding, AspectRatio, PlotStyle, Axes, AxesLabel, AxesOrigin, AxesStyle, TicksStyle, LabelStyle, Frame, FrameLabel, FrameStyle, FrameTicks, RotateLabel, GridLines, GridLinesStyle, Prolog, Epilog, PlotLabel, Background, ImageSize, ColorFunction (a function, or named ramp: "Rainbow"/"CoolTones"/"WarmTones"/"Greyscale"/"Temperature"), ColorFunctionScaling (default True), Filling (Axis/Bottom/Top/a number), FillingStyle, PlotLegends (Automatic/"Expressions"/an explicit list), RegionFunction, Exclusions.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Plot[Sin[x], {x, 0, 2 Pi}]
Out[1]= -Graphics-

In[2]:= Plot[1/x, {x, -2, 2}]
Out[2]= -Graphics-

In[3]:= Plot[Sin[x], {x, a, b}]
Out[3]= Plot[Sin[x], {x, a, b}]

In[4]:= Plot[{Sin[1/x], Cos[1/x]}, {x, -Pi, Pi}]
Out[4]= -Graphics-
```

### Options (4)

```mathematica
In[5]:= Plot[Sin[x] + Sin[7 x], {x, -2, 2}, Mesh -> All]
Out[5]= -Graphics-

In[6]:= Plot[Sin[x], {x, 0, 2 Pi}, GridLines -> Automatic, Epilog -> {Red, Point[{0, 0}]}, AxesOrigin -> {0, 0}]
Out[6]= -Graphics-

In[7]:= Plot[Sin[x], {x, 0, 2 Pi}, ColorFunction -> "Rainbow", Filling -> Axis]
Out[7]= -Graphics-

In[8]:= Plot[{Sin[x], Cos[x]}, {x, 0, 2 Pi}, PlotLegends -> "Expressions"]
Out[8]= -Graphics-
```

## Algorithm

plot.c — Plot[f, {x, xmin, xmax}, opts...].

HoldAll, like Table/Do: f and the iterator spec must not be pre-evaluated (x has no value yet). Splits trailing options into the sampler's own (PlotPoints/MaxRecursion/MaxPlotPoints/Mesh/RegionFunction/ Exclusions/ColorFunction/Filling/..., consumed here) and everything else (PlotRange/AspectRatio/PlotStyle/Axes/.../ImageSize, copied through onto the resulting Graphics[...] unevaluated -- render.c is the single place that interprets those, whether reached via Plot's auto-display or a later Show[]).

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## See also

[Show](../../graphics/Show/), [HoldAll](../../expression-information/HoldAll/)

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
- Tests: [`tests/test_graphics_sampling.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics_sampling.c)
