# ComplexPlot3D

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
ComplexPlot3D[f, {z, zmin, zmax}, opts...]
    Three-dimensional surface plot of a complex function: height = |f(z)|,
    colour = Arg(f(z)) via the thermal ramp (same default as Plot3D).
    z is bound to Complex[x, y] at each grid point; the result is a
    Graphics3D[...] object rendered in an orbit-camera window.
    ComplexPlot3D is HoldAll.

    Options:
      PlotPoints          grid resolution per axis (default 200)
      ColorFunction       f[re, im] → color, or named ramp string;
                          "PhaseRings" recommended (see ComplexPlot)
      ColorFunctionScaling True (default): scale re/im to [0,1]
      RegionFunction      f[x,y] mask
      PlotLegends         Automatic / True: attach a vertical phase color
                          scale bar (-π at bottom, π at top)
      Lighting -> None    disables Lambertian shading (recommended for
                          accurate phase colours; default Automatic)
      Standard Graphics3D options pass through to the result.

    Examples:
      ComplexPlot3D[z^2, {z, -2-2I, 2+2I}]
      ComplexPlot3D[Sin[z], {z, -2-2I, 2+2I}, Lighting->None]
      ComplexPlot3D[1/z, {z, -2-2I, 2+2I}, PlotLegends->Automatic, Lighting->None]
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
