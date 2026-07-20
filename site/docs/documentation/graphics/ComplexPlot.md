# ComplexPlot

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
ComplexPlot[f, {z, zmin, zmax}, opts...]
    Domain-colouring plot of the complex function f over the rectangular
    region in the complex plane with corners zmin and zmax.  z is bound
    to Complex[x, y] at each grid point; f must return a complex or real
    number.  The color of each cell encodes Arg(f(z)) via the thermal
    ramp (same default as DensityPlot), and brightness encodes |f(z)|.
    ComplexPlot is HoldAll: f and the iterator spec are held unevaluated
    until z is given a numeric complex value.

    Options:
      PlotPoints          grid resolution per axis (default 200)
      ColorFunction       f[re, im] → color, or a named ramp string.
                          Named ramps: "PhaseRings" (hue=phase,
                          brightness=log|w| rings — highlights poles
                          and zeros), "Rainbow", "CoolTones",
                          "WarmTones", "Greyscale", "Temperature"
                          (others keyed to normalised Arg)
      ColorFunctionScaling True (default): scale re/im to [0,1] before
                           calling a custom ColorFunction
      RegionFunction      f[x,y] mask; excluded cells are not drawn
      PlotLegends         Automatic / True: attach a vertical phase color
                          scale bar (thermal ramp, -π at bottom, π at top)
      Standard Graphics options (Axes, AspectRatio→1, Frame, PlotRange,
      AxesLabel, GridLines, ImageSize, Background, PlotLabel, …) pass
      through to the Graphics[...] result.

    Examples:
      ComplexPlot[z^2, {z, -2-2I, 2+2I}]
      ComplexPlot[Sin[z], {z, -Pi-Pi*I, Pi+Pi*I}]
      ComplexPlot[1/(z^2+1), {z, -2-2I, 2+2I}, PlotPoints->80]
      ComplexPlot[(z^2+1)/(z^2-1), {z, -2-2I, 2+2I}, PlotLegends->Automatic]
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
