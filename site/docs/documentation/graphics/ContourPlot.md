# ContourPlot

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
ContourPlot[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
    Generates iso-contour lines of f(x,y) using the marching squares algorithm
    and returns a Graphics[...] object (auto-displayed). ContourPlot is HoldAll:
    f is held unevaluated until x and y are bound to numeric values.

    Options:
      Contours         - Integer n (n evenly spaced auto levels, default 10), or
                         {c1, c2, ...} (explicit contour values).
      ContourStyle     - Style directive(s) for the contour lines. A single
                         directive is applied to all levels; a List cycles through
                         the levels. Automatic (default) colours by height.
                         None/False suppresses lines (leaves only shading).
      ContourLabels    - True: draw the z value at the midpoint of each level's
                         first visible segment. Default False.
      ContourShading   - True: fill each grid cell by its z value (via
                         ColorFunction or the built-in thermal gradient).
                         False/None: lines only. Automatic (default): shade when
                         ColorFunction is set, otherwise lines only.
      ColorFunction    - A function f[t] → color (t in [0,1] after scaling), or
                         a string: "Rainbow" (Hue ramp) or "Temperature" (blue-
                         cyan-yellow-red). Applied to shading and auto line colors.
      ColorFunctionScaling - True (default): normalise z to [0,1] before calling
                             ColorFunction. False: pass raw z.
      PlotPoints       - Grid resolution per axis (default 25; increase for
                         smoother contours).
      RegionFunction   - f[x,y] mask: cells where the function is False are
                         skipped (neither shaded nor contoured).
      Standard Graphics options (Axes, AspectRatio, Frame, PlotRange,
      AxesLabel, GridLines, ImageSize, Background, PlotLabel, Prolog,
      Epilog, ...) pass through to the Graphics[...] result.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

- `HoldAll`, `Protected`.
- Declines to evaluate if bounds aren't numeric or the function argument

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
