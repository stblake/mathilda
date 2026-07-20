# Plot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Plot[f, {x, xmin, xmax}, opts...]
    Adaptively samples f over [xmin, xmax], displays the resulting curve in an interactive window, and returns it as a Graphics[...] object. A list of functions Plot[{f1, f2, ...}, {x, xmin, xmax}] draws each on the same axes in a distinct palette colour. Options: PlotPoints (initial sample count, default 50), MaxRecursion (adaptive refinement depth, default 6), MaxPlotPoints (overall point cap, default Infinity), Mesh (All overlays the evaluation points as dots; default None), PlotRange, PlotRangePadding, AspectRatio, PlotStyle, Axes, AxesLabel, AxesOrigin, AxesStyle, TicksStyle, LabelStyle, Frame, FrameLabel, FrameStyle, FrameTicks, RotateLabel, GridLines, GridLinesStyle, Prolog, Epilog, PlotLabel, Background, ImageSize, ColorFunction (a function, or named ramp: "Rainbow"/"CoolTones"/"WarmTones"/"Greyscale"/"Temperature"), ColorFunctionScaling (default True), Filling (Axis/Bottom/Top/a number), FillingStyle, PlotLegends (Automatic/"Expressions"/an explicit list), RegionFunction, Exclusions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Plot[Sin[x], {x, a, b}]
Out[1]= Plot[Sin[x], {x, a, b}]
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
