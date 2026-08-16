# ArrayPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ArrayPlot[array, opts...]`**

Renders a 2D array (nested List or NDArray) as a grid of coloured cells, one per array entry -- a discrete heatmap with no interpolation between cells. Row 1 of array is drawn at the top, column 1 at the left. Not HoldAll: array is an ordinary evaluated expression. Returns a Graphics\[...\] object (auto-displayed). ArrayPlot\[colorArray, opts...\] If a cell is already a colour literal (RGBColor/GrayLevel/Hue/ CMYKColor), it paints that colour directly instead of one derived from ColorFunction -- ArrayPlot doubles as a raw pixel-grid renderer, e.g. ArrayPlot\[{{Red, Blue}, {Blue, Red}}\]. Numeric and colour cells freely mix within the same array: ArrayPlot\[{{1, 0, Pink}, {0, 1, Red}}\] calls out two cells explicitly while the rest still follow the normal heatmap. Options: ColorFunction        named ramp string or f\[t\]-\>color (t in \[0,1\]). Ramps: "Greyscale" (default: white low, black high -- matches Mathematica's ArrayPlot), "Rainbow", "Temperature", "CoolTones", "WarmTones", all keyed to the normalised entry value. Only applies to numeric cells. ColorFunctionScaling True (default): normalise entries to \[0,1\] before calling ColorFunction; False: raw value ColorRules           {v1 -\> c1, v2 -\> c2, ...} (or a single v -\> c): an explicit colour for numeric cells whose value exactly equals v, checked before ColorFunction. Cells matching no rule still get the normal scaled ColorFunction colour. Mesh                 All/True: draw grey grid lines between cells; None (default): no lines PlotLegends          Automatic: attach a vertical colour scale bar (only when at least one cell is numeric) Standard Graphics options (Axes, AspectRatio -\> rows/cols by default, Frame, PlotRange, ImageSize, Background, PlotLabel, ...) pass through to the Graphics\[...\] result. Examples: ArrayPlot\[{{1, 0, 1}, {0, 1, 0}, {1, 0, 1}}\] ArrayPlot\[RandomReal\[1, {20, 20}\], ColorFunction -\> "Rainbow"\] ArrayPlot\[Table\[Mod\[i + j, 2\], {i, 10}, {j, 10}\], Mesh -\> All\] ArrayPlot\[{{1, 0}, {0, 1}}, ColorRules -\> {1 -\> Pink, 0 -\> Yellow}\] ArrayPlot\[{{1, 0, 0, Pink}, {1, 1, 0, Pink}, {1, 0, 1, Red}}\]

## Examples (6)

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

## Algorithm

arrayplot.c — ArrayPlot[array, opts...]

Renders a 2D array of values as a grid of coloured cells (a discrete heatmap, no interpolation): row 0 of the array is drawn at the top, column 0 at the left, matching Mathematica's ArrayPlot orientation. Each numeric cell's value is normalised to [0,1] (ColorFunctionScaling, default True) and mapped through ColorFunction. The default ramp is the shared "Greyscale" ramp from plot_common (white at the minimum, black at the maximum) -- the same ramp DensityPlot/ComplexPlot expose by name, and Mathematica's own ArrayPlot default.

A cell that is already a colour literal (RGBColor/GrayLevel/Hue/ CMYKColor) is painted directly instead, and numeric and colour cells freely mix within the same array -- ArrayPlot[{{1, 0, Pink}, {0, 1, Red}}] calls out two cells explicitly while the rest still follow the normal heatmap, and ArrayPlot[{{Red, Blue}, {Blue, Red}}] (every cell a colour) lets ArrayPlot double as a raw pixel-grid renderer.

Unlike the HoldAll function-sampling plots (DensityPlot, Plot, ...), ArrayPlot is not HoldAll: its argument is an ordinary already-evaluated nested List or NDArray. An NDArray input (a dense numeric buffer by construction) goes through na_load_matrix (src/linalg/numarray.h); a nested List goes through arrayplot_load below, which classifies each cell individually.

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [NDArray](../../linear-algebra/NDArray/), [DensityPlot](../../graphics/DensityPlot/), [HoldAll](../../expression-information/HoldAll/), [CMYKColor](../../graphics/CMYKColor/), [AspectRatio](../../other-advanced/AspectRatio/), [Frame](../../other-advanced/Frame/), [ImageSize](../../other-advanced/ImageSize/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
