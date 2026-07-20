# PolarPlot

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
PolarPlot[r, {theta, tmin, tmax}, opts...]
    Plots the polar curve r(theta) by converting to Cartesian
    coordinates {r*Cos[theta], r*Sin[theta]} and sampling adaptively
    over [tmin, tmax]. Returns a Graphics[...] object (auto-displayed).
PolarPlot[{r1, r2, ...}, {theta, tmin, tmax}, opts...]
    Multiple polar curves in distinct palette colours.

    Negative r values are plotted in the opposite direction (standard
    polar convention). Default AspectRatio -> 1 (equal axes).

    Options (same as ParametricPlot):
      PlotPoints          - initial sample count per curve (default 75)
      MaxRecursion        - adaptive refinement depth (default 6)
      MaxPlotPoints       - total point cap (default Infinity)
      Mesh                - All/True: overlay evaluation dots; None (default)
      ColorFunction       - f[t] or "Rainbow" (sweeps scaled theta)
      ColorFunctionScaling - True (default): normalise theta to [0,1]
      RegionFunction      - f[x,y] mask
      PlotStyle           - color/style directive(s)
      PlotLegends         - Automatic / "Expressions" / label list
      PolarAxes           - option keyword (accepted; polar grid overlay
                            not yet rendered)
      Standard Graphics options (AspectRatio, Axes, PlotRange,
      AxesLabel, Frame, GridLines, PlotLabel, Background, ImageSize,
      Prolog, Epilog) pass through to the Graphics[...] result.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
