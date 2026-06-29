# ListPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ListPlot[{y1, ..., yn}, opts...]
    Plots the values as points {i, yi} (a scatter/point plot). ListPlot[{{x1,y1}, ...}] plots the given coordinate pairs; ListPlot[{data1, data2, ...}] overlays each dataset in a distinct palette colour. Returns a Graphics[...] object. Options: Joined (connect points; default False), DataRange (x-range for heights), Filling (Axis/Bottom/Top/a number — draws stems), FillingStyle, PlotMarkers, PlotStyle, PlotLegends, and the Graphics options PlotRange, Axes (default True), AspectRatio (default 1/GoldenRatio), Frame, AxesLabel, GridLines, ImageSize, Background, PlotLabel, Prolog, Epilog.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
