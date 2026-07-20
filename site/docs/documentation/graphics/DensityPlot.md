# DensityPlot

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
DensityPlot[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
    Renders f(x,y) as a heatmap: each grid cell is coloured by its
    function value via ColorFunction (default: thermal blue→yellow ramp).
    DensityPlot is HoldAll: f is held unevaluated until x and y are
    bound to numeric values. Returns a Graphics[...] object.

    Options:
      PlotPoints          grid resolution per axis (default 50)
      ColorFunction       named ramp string or f[t]→color (t in [0,1]).
                          Ramps: "Rainbow", "CoolTones", "WarmTones",
                          "Greyscale", "Temperature" (all keyed to
                          normalised z value, t∈[0,1])
      ColorFunctionScaling True (default): normalise z to [0,1] before
                           calling ColorFunction; False: pass raw z
      RegionFunction      f[x,y] mask; excluded cells are not drawn
      PlotLegends         Automatic: attach a vertical color scale bar
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
