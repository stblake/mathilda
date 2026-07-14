# ParametricPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ParametricPlot[{fx, fy}, {t, tmin, tmax}, opts...]
    Adaptively samples the parametric curve (fx(t), fy(t)) over [tmin, tmax] and returns a Graphics[...] object (auto-displayed). The body may be any expression that evaluates to a 2-element {x,y} list (not just a literal {fx,fy}). Multiple curves: ParametricPlot[{{fx1,fy1}, ...}, {t,...}]. Two-iterator (filled region) form: ParametricPlot[body, {t,...}, {r,...}] samples a PlotPoints x PlotPoints grid and emits Polygon[] quads. Default AspectRatio -> 1 (both axes equally important). Options: PlotPoints (default 25), MaxRecursion (default 6), MaxPlotPoints, Mesh (All: dots for curves, grid lines for regions), PlotLegends (Automatic/"Expressions"/{labels...}: draws a legend), ColorFunction ("Rainbow" or f[t] / f[t,r]), ColorFunctionScaling (default True), RegionFunction (f[x,y] mask), PlotStyle, AspectRatio, Axes, PlotRange, PlotRangePadding, AxesLabel, AxesOrigin, Frame, FrameLabel, GridLines, Prolog, Epilog, PlotLabel, Background, ImageSize (all passed through to Graphics).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `HoldAll`, `Protected`.
- Declines to evaluate if bounds aren't numeric or `PlotPoints < 2`.
- Auto-displays exactly like `Graphics`/`Plot`.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
