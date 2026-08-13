# ListPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ListPlot[{y1, ..., yn}, opts...]`**

Plots the values as points {i, yi} (a scatter/point plot). ListPlot\[{{x1,y1}, ...}\] plots the given coordinate pairs; ListPlot\[{data1, data2, ...}\] overlays each dataset in a distinct palette colour. Returns a Graphics\[...\] object. Options: Joined (connect points; default False), DataRange (x-range for heights), Filling (Axis/Bottom/Top/a number — draws stems), FillingStyle, PlotMarkers, PlotStyle, PlotLegends, and the Graphics options PlotRange, Axes (default True), AspectRatio (default 1/GoldenRatio), Frame, AxesLabel, GridLines, ImageSize, Background, PlotLabel, Prolog, Epilog.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= ListPlot[{1, 4, 9, 16, 25}]
Out[1]= -Graphics-

In[2]:= ListPlot[{{0, 0}, {1, 1}, {2, 4}, {3, 9}}]
Out[2]= -Graphics-
```

### Options (3)

```mathematica
In[3]:= ListPlot[{{1, 1}, {2, 4}, {3, 9}}, Joined -> True]
Out[3]= -Graphics-

In[4]:= ListPlot[{Table[Sin[n], {n, 20}], Table[Cos[n], {n, 20}]}, PlotLegends -> {"sin", "cos"}]
Out[4]= -Graphics-

In[5]:= ListPlot[{1, 4, 9, 16}, Filling -> Axis]
Out[5]= -Graphics-
```

## Algorithm

listplot.c — ListPlot[data, opts...].

Unlike Plot (which is HoldAll because its function body must stay symbolic while x is unbound), ListPlot's data is concrete and must be evaluated (so ListPlot[Table[i^2, {i, 5}]] / ListPlot[Range[10]] work), so ListPlot is a plain protected builtin. Its arguments — the data and the option values — therefore arrive already evaluated (named colours like Red are RGBColor[...] by the time we see them, exactly as a bare Graphics[]'s own arguments would be).

The work is purely constructive: classify the data into one or more datasets, turn each into Point[...] (or Line[...] when Joined) primitives, and wrap them in a Graphics[...] carrying the passthrough options. The existing renderer (render.c) interprets PlotRange/Axes/AspectRatio/Frame/ PlotStyle/PlotLegends, so ListPlot inherits all of them for free. Sampler helpers numericize_bound/palette_color are shared from plot.c (see plot.h).

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Plot](../../graphics/Plot/), [Show](../../graphics/Show/), [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
