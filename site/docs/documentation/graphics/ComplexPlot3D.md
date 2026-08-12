# ComplexPlot3D

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`ComplexPlot3D[f, {z, zmin, zmax}, opts...]`**

Three-dimensional surface plot of a complex function: height = |f(z)|, colour = Arg(f(z)) via the thermal ramp (same default as Plot3D). z is bound to Complex\[x, y\] at each grid point; the result is a Graphics3D\[...\] object rendered in an orbit-camera window. ComplexPlot3D is HoldAll. Options: PlotPoints          grid resolution per axis (default 200) ColorFunction       f\[re, im\] → color, or named ramp string; "PhaseRings" recommended (see ComplexPlot) ColorFunctionScaling True (default): scale re/im to \[0,1\] RegionFunction      f\[x,y\] mask PlotLegends         Automatic / True: attach a vertical phase color scale bar (-π at bottom, π at top) Lighting -\> None    disables Lambertian shading (recommended for accurate phase colours; default Automatic) Standard Graphics3D options pass through to the result. Examples: ComplexPlot3D\[z^2, {z, -2-2I, 2+2I}\] ComplexPlot3D\[Sin\[z\], {z, -2-2I, 2+2I}, Lighting-\>None\] ComplexPlot3D\[1/z, {z, -2-2I, 2+2I}, PlotLegends-\>Automatic, Lighting-\>None\]

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

## See also

[Plot3D](../../graphics/Plot3D/), [HoldAll](../../expression-information/HoldAll/), [ComplexPlot](../../graphics/ComplexPlot/), [Lighting](../../other-advanced/Lighting/)

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
