# ComplexPlot3D

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`ComplexPlot3D[f, {z, zmin, zmax}, opts...]`**

Three-dimensional surface plot of a complex function: height = |f(z)|, colour = Arg(f(z)) via the thermal ramp (same default as Plot3D). z is bound to Complex\[x, y\] at each grid point; the result is a Graphics3D\[...\] object rendered in an orbit-camera window. ComplexPlot3D is HoldAll. Options: PlotPoints          grid resolution per axis (default 200) ColorFunction       a named ramp string, or a function of the eight arguments Re\[z\], Im\[z\], Abs\[z\], Arg\[z\], Re\[f\], Im\[f\], Abs\[f\], Arg\[f\] (see ComplexPlot); "PhaseRings" recommended ColorFunctionScaling True (default): scale the eight args to \[0,1\] RegionFunction      f\[x,y\] mask PlotLegends         Automatic / True: attach a vertical phase color scale bar (-π at bottom, π at top) Lighting -\> None    disables Lambertian shading (recommended for accurate phase colours; default Automatic) Standard Graphics3D options pass through to the result. Examples: ComplexPlot3D\[z^2, {z, -2-2I, 2+2I}\] ComplexPlot3D\[Sin\[z\], {z, -2-2I, 2+2I}, Lighting-\>None\] ComplexPlot3D\[1/z, {z, -2-2I, 2+2I}, PlotLegends-\>Automatic, Lighting-\>None\]

## Examples

_No verified examples yet for this function._

## Algorithm

complexplot.c — ComplexPlot and ComplexPlot3D.

Both functions share the same domain-parsing logic and option set. The only structural difference is what they build from the evaluated grid: ComplexPlot emits Rectangle primitives into Graphics[], and ComplexPlot3D emits Polygon quads (height = |w|, colour = arg(w)) into Graphics3D[].

Coloring convention: the default maps the phase arg(w) onto the cyclic "Cyclic" ramp (t = (atan2(im, re) + π) / (2π) ∈ [0, 1]) and folds the modulus in as an HSL lightness — zeros fade to black, poles to white — so ComplexPlot[f] and ComplexPlot[f, ColorFunction -> "Cyclic"] are identical. A custom ColorFunction receives the eight Mathematica arguments

```text
  Re[z], Im[z], Abs[z], Arg[z], Re[f], Im[f], Abs[f], Arg[f]
```

(so #8 is the phase of the value); with ColorFunctionScaling→True (default) each is scaled to [0,1] across the sampled domain.

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

**See also:** [Plot3D](../../graphics/Plot3D/), [HoldAll](../../expression-information/HoldAll/), [ComplexPlot](../../graphics/ComplexPlot/), [Lighting](../../other-advanced/Lighting/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
