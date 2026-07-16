# StreamPlot

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
StreamPlot[{vx, vy}, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
    Traces streamlines of the 2-D vector field {vx, vy} by RK4 integration
    from a grid of seed points, and returns a Graphics[{Arrow[...], ...}, opts]
    object (auto-displayed). StreamPlot is HoldAll: vx, vy, and the iterator
    specs are held unevaluated until x and y are given numeric values.

    Options:
      StreamPoints  - Integer n (n x n seed grid) or Automatic (default 15 x 15).
      StreamScale   – Automatic (8%% of domain diagonal, default), None (full run),
                       or a real fraction of the domain diagonal.
      StreamStyle   – Style directive(s) applied to all streams.
      StreamColorFunction / ColorFunction
                     – f[x,y,vx,vy,speed] (or fewer args) returning a color,
                       or "Rainbow" (hue = scaled speed).
      RegionFunction – f[x,y] mask; seeds outside the region are skipped.
      PlotLegends   – Automatic / "Expressions" / explicit label list.
      Standard Graphics options (PlotRange, Axes, AspectRatio, Frame, …)
                       pass through to the Graphics[...] result.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `HoldAll`, `Protected`.
- RK4 integration; step size adapts to seed density and domain size.
- Declines to evaluate if the field arg is not a 2-element List, or if

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
