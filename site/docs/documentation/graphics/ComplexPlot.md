# ComplexPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ComplexPlot[f, {z, zmin, zmax}, opts...]`**

Domain-colouring plot of the complex function f over the rectangular region in the complex plane with corners zmin and zmax.  z is bound to Complex\[x, y\] at each grid point; f must return a complex or real number.  The color of each cell encodes Arg(f(z)) via the thermal ramp (same default as DensityPlot), and brightness encodes |f(z)|. ComplexPlot is HoldAll: f and the iterator spec are held unevaluated until z is given a numeric complex value. Options: PlotPoints          grid resolution per axis (default 200) ColorFunction       f\[re, im\] → color, or a named ramp string. Named ramps: "PhaseRings" (hue=phase, brightness=log|w| rings — highlights poles and zeros), "Rainbow", "CoolTones", "WarmTones", "Greyscale", "Temperature" (others keyed to normalised Arg) ColorFunctionScaling True (default): scale re/im to \[0,1\] before calling a custom ColorFunction RegionFunction      f\[x,y\] mask; excluded cells are not drawn PlotLegends         Automatic / True: attach a vertical phase color scale bar (thermal ramp, -π at bottom, π at top) Standard Graphics options (Axes, AspectRatio→1, Frame, PlotRange, AxesLabel, GridLines, ImageSize, Background, PlotLabel, …) pass through to the Graphics\[...\] result. Examples: ComplexPlot\[z^2, {z, -2-2I, 2+2I}\] ComplexPlot\[Sin\[z\], {z, -Pi-Pi\*I, Pi+Pi\*I}\] ComplexPlot\[1/(z^2+1), {z, -2-2I, 2+2I}, PlotPoints-\>80\] ComplexPlot\[(z^2+1)/(z^2-1), {z, -2-2I, 2+2I}, PlotLegends-\>Automatic\]

## Examples

_No verified examples yet for this function._

## Algorithm

complexplot.c — ComplexPlot and ComplexPlot3D.

Both functions share the same domain-parsing logic and option set. The only structural difference is what they build from the evaluated grid: ComplexPlot emits Rectangle primitives into Graphics[], and ComplexPlot3D emits Polygon quads (height = |w|, colour = arg(w)) into Graphics3D[].

Coloring convention: the default phase-to-color mapping uses the same thermal_rgb ramp that DensityPlot and Plot3D use, keyed to the normalised argument: t = (atan2(im, re) + π) / (2π) ∈ [0, 1]. This keeps the palette consistent across the whole graphics engine. A custom ColorFunction receives (re, im) as a 2-arg call; with ColorFunctionScaling→True (default) the arguments are first scaled so that re ∈ [0,1] and im ∈ [0,1] across the sampled domain.

Both are HoldAll: the body and the iterator spec are held unevaluated until z is bound to Complex[x, y] at each grid point.

Domain spec:

```text
  {z, zmin, zmax}  — z is the complex iterator variable; zmin and zmax
  are evaluated and their Re/Im parts define the rectangular plotting
  domain: xmin=Re(zmin), xmax=Re(zmax), ymin=Im(zmin), ymax=Im(zmax).
  Both endpoints may be real (imaginary part = 0). 
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [DensityPlot](../../graphics/DensityPlot/), [HoldAll](../../expression-information/HoldAll/), [Arg](../../arithmetic/Arg/), [AspectRatio](../../other-advanced/AspectRatio/), [Frame](../../other-advanced/Frame/), [ImageSize](../../other-advanced/ImageSize/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
