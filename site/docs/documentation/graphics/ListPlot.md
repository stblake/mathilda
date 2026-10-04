# ListPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ListPlot[{y1, ..., yn}, opts...]`**

Plots the values as points {i, yi} (a scatter/point plot). ListPlot\[{{x1,y1}, ...}\] plots the given coordinate pairs; ListPlot\[{data1, data2, ...}\] overlays each dataset in a distinct palette colour. Returns a Graphics\[...\] object. Options: Joined (connect points; default False), DataRange (x-range for heights), Filling (Axis/Bottom/Top/a number — draws stems), FillingStyle, PlotMarkers, PlotStyle, PlotLegends, and the Graphics options PlotRange, Axes (default True), AspectRatio (default 1/GoldenRatio), Frame, AxesLabel, GridLines, ImageSize, Background, PlotLabel, Prolog, Epilog.

## Examples (10)

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

### Applications (5)

```mathematica
In[6]:= ListPlot[{1, 4, 9, 16, 25}]
Out[6]= -Graphics-

In[7]:= ListPlot[{{0, 0}, {1, 1}, {2, 4}, {3, 9}}]
Out[7]= -Graphics-

In[8]:= ListPlot[Table[Sin[n], {n, 20}], Joined -> True]
Out[8]= -Graphics-

In[9]:= Length[Cases[ListPlot[Range[10]], _Point, Infinity]]
Out[9]= 1

In[10]:= Options[ListPlot, AspectRatio]
Out[10]= {AspectRatio -> 1/GoldenRatio}
```

## Algorithm

listplot.c — ListPlot[data, opts...].

Unlike Plot (which is HoldAll because its function body must stay symbolic while x is unbound), ListPlot's data is concrete and must be evaluated (so ListPlot[Table[i^2, {i, 5}]] / ListPlot[Range[10]] work), so ListPlot is a plain protected builtin. Its arguments — the data and the option values — therefore arrive already evaluated (named colours like Red are RGBColor[...] by the time we see them, exactly as a bare Graphics[]'s own arguments would be).

The work is purely constructive: classify the data into one or more datasets, turn each into Point[...] (or Line[...] when Joined) primitives, and wrap them in a Graphics[...] carrying the passthrough options. The existing renderer (render.c) interprets PlotRange/Axes/AspectRatio/Frame/ PlotStyle/PlotLegends, so ListPlot inherits all of them for free. Sampler helpers numericize_bound/palette_color are shared from plot.c (see plot.h).

## Implementation notes

**Algorithm.** `builtin_listplot` is a plain `Protected` builtin — **not**
`HoldAll` — so its data is already evaluated (`ListPlot[Range[10]]`,
`ListPlot[Table[i^2, {i, 5}]]` and named colours in options all work). `classify`
/ `extract_dataset` sort the argument into one or more datasets: a flat list with
no sublists is read as **heights** `{y_i}` plotted at `x` running over `DataRange`
(default `1..n`), with non-numeric entries skipped but still consuming their index
slot; a list in which *every* element is a numeric 2-pair is read as **explicit
points** `{{x, y}, ...}`; anything else (some elements are sublists) is read as
**multiple datasets**, one per sublist, in distinct `palette_color`s.
`DataRange -> All` forces the multi reading of a pair list. Numeric coercion goes
through `to_double` (literal `Integer`/`Real` fast path, else `N[]`). Each dataset
is emitted by `emit_dataset` as a coordinate `List` wrapped in `Line[...]` when
`Joined -> True` (default `False`) and otherwise `Point[...]`; multi-dataset
styling gets a per-dataset palette colour with non-colour directives in their own
`List` scope. `Filling` draws continuous quads under a joined curve
(`gfx_build_fill_quads`) or vertical stems under unjoined points;
`ScalingFunctions` transform world coordinates in place; `PlotLegends` attaches
legend metadata. The result is an inert `Graphics[prims, opts...]`.

**Data structures.** `PointSet {double* xs; double* ys; size_t n}`, one per live
dataset (`classify` returns an array of borrowed dataset `Expr*`); the option
bundle from `split_options`.

**Complexity / limits.** `O(total points)` — there is no sampling, the data is
plotted as given. Default `Axes -> True`, `AspectRatio -> 1/GoldenRatio`
(distinct from the field plotters' `1`), `PlotStyle` colour
`RGBColor[0.2, 0.4, 0.8]`. The default marker size is resolved by the renderer
from scatter density rather than emitted here; `PlotMarkers` is accepted but its
glyphs are not yet drawn.

**Attributes:** `Protected`.

## References

**See also:** [Plot](../../graphics/Plot/), [Show](../../graphics/Show/), [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/graphics/listplot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/listplot.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)

## Notes & additional examples

### Notes

`ListPlot` is `Protected` but **not** `HoldAll`, so its data is evaluated before
classification — `ListPlot[Range[10]]` and `ListPlot[Table[i^2, {i, 5}]]` both
work. The argument is read as heights (a flat list, plotted at `{i, y_i}` over
`DataRange`), explicit coordinate pairs (every element a numeric 2-list), or
several datasets (some elements are sublists), each in a distinct `ColorData[97]`
palette colour. Non-numeric entries are skipped but still consume their index
slot.

`Joined -> True` emits one `Line` per dataset instead of a `Point` cloud;
`Filling` draws continuous quads under a joined curve or vertical stems under
points. The default `AspectRatio` is `1/GoldenRatio` (the field plotters default
to `1`), and the point marker size is chosen by the renderer from scatter density
unless an explicit `PointSize` overrides it.
