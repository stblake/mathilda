# ArrayPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ArrayPlot[array, opts...]`**

Renders a 2D array (nested List or NDArray) as a grid of coloured cells, one per array entry -- a discrete heatmap with no interpolation between cells. Row 1 of array is drawn at the top, column 1 at the left. Not HoldAll: array is an ordinary evaluated expression. Returns a Graphics\[...\] object (auto-displayed). ArrayPlot\[colorArray, opts...\] If a cell is already a colour literal (RGBColor/GrayLevel/Hue/ CMYKColor), it paints that colour directly instead of one derived from ColorFunction -- ArrayPlot doubles as a raw pixel-grid renderer, e.g. ArrayPlot\[{{Red, Blue}, {Blue, Red}}\]. Numeric and colour cells freely mix within the same array: ArrayPlot\[{{1, 0, Pink}, {0, 1, Red}}\] calls out two cells explicitly while the rest still follow the normal heatmap. Options: ColorFunction        named ramp string or f\[t\]-\>color (t in \[0,1\]). Ramps: "Greyscale" (default: white low, black high -- matches Mathematica's ArrayPlot), "Rainbow", "Temperature", "CoolTones", "WarmTones", all keyed to the normalised entry value. Only applies to numeric cells. ColorFunctionScaling True (default): normalise entries to \[0,1\] before calling ColorFunction; False: raw value ColorRules           {v1 -\> c1, v2 -\> c2, ...} (or a single v -\> c): an explicit colour for numeric cells whose value exactly equals v, checked before ColorFunction. Cells matching no rule still get the normal scaled ColorFunction colour. Mesh                 All/True: draw grey grid lines between cells; None (default): no lines PlotLegends          Automatic: attach a vertical colour scale bar (only when at least one cell is numeric) Standard Graphics options (Axes, AspectRatio -\> rows/cols by default, Frame, PlotRange, ImageSize, Background, PlotLabel, ...) pass through to the Graphics\[...\] result. Examples: ArrayPlot\[{{1, 0, 1}, {0, 1, 0}, {1, 0, 1}}\] ArrayPlot\[RandomReal\[1, {20, 20}\], ColorFunction -\> "Rainbow"\] ArrayPlot\[Table\[Mod\[i + j, 2\], {i, 10}, {j, 10}\], Mesh -\> All\] ArrayPlot\[{{1, 0}, {0, 1}}, ColorRules -\> {1 -\> Pink, 0 -\> Yellow}\] ArrayPlot\[{{1, 0, 0, Pink}, {1, 1, 0, Pink}, {1, 0, 1, Red}}\]

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ArrayPlot[{{1, 0, 1}, {0, 1, 0}, {1, 0, 1}}]
Out[1]= -Graphics-

In[2]:= ArrayPlot[{{Red, Blue}, {Blue, Red}}]
Out[2]= -Graphics-

In[3]:= ArrayPlot[{{1, 0, 0, Pink}, {1, 1, 0, Pink}, {1, 0, 1, Red}}]
Out[3]= -Graphics-
```

### Options (3)

```mathematica
In[4]:= ArrayPlot[RandomReal[1, {20, 20}], ColorFunction -> "Rainbow", Mesh -> All]
Out[4]= -Graphics-

In[5]:= ArrayPlot[Table[Mod[i + j, 5], {i, 10}, {j, 10}], PlotLegends -> Automatic]
Out[5]= -Graphics-

In[6]:= ArrayPlot[{{1, 0, 0.5}, {0, 1, 0.5}}, ColorRules -> {1 -> Pink, 0 -> Yellow}]
Out[6]= -Graphics-
```

### Applications (5)

```mathematica
In[7]:= ArrayPlot[{{1, 0, 1}, {0, 1, 0}, {1, 0, 1}}]
Out[7]= -Graphics-

In[8]:= ArrayPlot[{{Red, Blue}, {Blue, Red}}]
Out[8]= -Graphics-

In[9]:= ArrayPlot[{{1, 0, 0.5}, {0, 1, 0.5}}, ColorRules -> {1 -> Pink, 0 -> Yellow}]
Out[9]= -Graphics-

In[10]:= Length[Cases[ArrayPlot[{{1, 0}, {0, 1}}], _Rectangle, Infinity]]
Out[10]= 4

In[11]:= Attributes[ArrayPlot]
Out[11]= {Protected}
```

## Algorithm

arrayplot.c — ArrayPlot[array, opts...]

Renders a 2D array of values as a grid of coloured cells (a discrete heatmap, no interpolation): row 0 of the array is drawn at the top, column 0 at the left, matching Mathematica's ArrayPlot orientation. Each numeric cell's value is normalised to [0,1] (ColorFunctionScaling, default True) and mapped through ColorFunction. The default ramp is the shared "Greyscale" ramp from plot_common (white at the minimum, black at the maximum) -- the same ramp DensityPlot/ComplexPlot expose by name, and Mathematica's own ArrayPlot default.

A cell that is already a colour literal (RGBColor/GrayLevel/Hue/ CMYKColor) is painted directly instead, and numeric and colour cells freely mix within the same array -- ArrayPlot[{{1, 0, Pink}, {0, 1, Red}}] calls out two cells explicitly while the rest still follow the normal heatmap, and ArrayPlot[{{Red, Blue}, {Blue, Red}}] (every cell a colour) lets ArrayPlot double as a raw pixel-grid renderer.

Unlike the HoldAll function-sampling plots (DensityPlot, Plot, ...), ArrayPlot is not HoldAll: its argument is an ordinary already-evaluated nested List or NDArray. An NDArray input (a dense numeric buffer by construction) goes through na_load_matrix (src/linalg/numarray.h); a nested List goes through arrayplot_load below, which classifies each cell individually.

## Implementation notes

**Algorithm.** `builtin_arrayplot` is a plain `Protected` builtin — **not**
`HoldAll` — so its argument is an already-evaluated array, rendered as a discrete
grid of coloured cells with **no interpolation** (unlike `DensityPlot`). Two
ingestion paths: an `EXPR_NDArray` is loaded by `na_load_matrix`
(`src/linalg/numarray.h`, a dense numeric buffer, every cell numeric); a nested
`List` goes to `arrayplot_load`, which walks the rectangular list-of-lists and
classifies each cell as either a **colour literal**
(`RGBColor`/`GrayLevel`/`Hue`/`CMYKColor`, painted directly) or a **plain number**
(`expr_to_real_double`). Numeric and colour cells may mix freely; a ragged shape,
a non-list-of-lists, or a cell that is neither colour nor number makes the whole
call decline (`NULL`). Each numeric cell is coloured by `ap_color`: the default
`ColorFunction` is **"Greyscale"** (white at the array minimum, black at the
maximum — matching Mathematica), with `ColorFunctionScaling -> True` normalising
`(v - zmin)/zspan`; `ColorRules -> {key -> colour}` overrides exact-matching
values before the ramp. Row 0 is drawn at the **top** (`y0 = rows-1-i`), column 0
at the **left**; each cell is a colour directive plus a `Rectangle`, overlapping
its later-drawn neighbours to close seams. `Mesh -> All` adds `rows+1` horizontal
and `cols+1` vertical `GrayLevel[0.5]` grid `Line`s. `PlotLegends -> Automatic`
appends a `$StreamColorBar` only when at least one cell is numeric. The result is
an inert `Graphics[prims, opts...]`.

**Data structures.** `ArrayPlotOpts`; parallel arrays `const Expr** cell_colors`
(borrowed colour-literal pointers, or NULL) and `double* buf` (values), each
`rows*cols`; `Expr** prims` sized `rows*cols*2 + mesh + 4`. `AspectRatio`
defaults to `rows/cols` (square cells); `PlotRange` is the exact extent
`{{0, cols}, {0, rows}}`.

**Complexity / limits.** `O(rows*cols)`, no sampling — a literal heatmap or, when
every cell is a colour, a raw pixel-grid renderer. Defaults `Frame -> True`,
`Axes -> False`.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [NDArray](../../linear-algebra/NDArray/), [DensityPlot](../../graphics/DensityPlot/), [HoldAll](../../expression-information/HoldAll/), [CMYKColor](../../graphics/CMYKColor/), [AspectRatio](../../other-advanced/AspectRatio/), [Frame](../../other-advanced/Frame/), [ImageSize](../../other-advanced/ImageSize/)

- Source: [`src/graphics/arrayplot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/arrayplot.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)

## Notes & additional examples

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
