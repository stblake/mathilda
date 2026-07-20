# VectorPlot

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
VectorPlot[{vx, vy}, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
    Draws a grid of arrows showing the direction (and optionally
    magnitude) of the vector field {vx, vy} at each grid point.
    VectorPlot is HoldAll: vx, vy are held unevaluated until x and y
    are bound to numeric values. Returns a Graphics[...] object.

    Options:
      VectorPoints   integer n → n×n grid (default 15); Automatic = 15
      VectorScale    Automatic: equal-length arrows (direction only)
                     None: proportional to magnitude
                     real f: arrow length = f × grid spacing
      VectorStyle    style directive(s) applied to all arrows
      ColorFunction  named ramp string (keyed to speed) or
                     f[vx,vy,speed]/f[speed]→color.
                     Ramps: "Rainbow", "CoolTones", "WarmTones",
                     "Greyscale", "Temperature"
      ColorFunctionScaling  True (default): normalise speed to [0,1]
      RegionFunction f[x,y] mask: skip grid points outside the region
      Standard Graphics options (Axes, AspectRatio→1, Frame, PlotRange,
      AxesLabel, GridLines, ImageSize, Background, PlotLabel, …) pass
      through to the Graphics[...] result.
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
