# ComplexPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ComplexPlot[f, {z, zmin, zmax}, opts...]`**

Domain-colouring plot of the complex function f over the rectangular region in the complex plane with corners zmin and zmax.  z is bound to Complex\[x, y\] at each grid point; f must return a complex or real number.  Each cell's hue encodes Arg(f(z)) on the cyclic "Cyclic" ramp, and |f(z)| sets an HSL lightness: zeros fade to black, poles to white.  ComplexPlot is HoldAll: f and the iterator spec are held unevaluated until z is given a numeric complex value. Options: PlotPoints          grid resolution per axis (default 400) ColorFunction       a named ramp string, or a function of the eight Mathematica arguments Re\[z\], Im\[z\], Abs\[z\], Arg\[z\], Re\[f\], Im\[f\], Abs\[f\], Arg\[f\] (so #8 is the phase of the value).  Named ramps: "PhaseRings" (hue=phase, brightness=log|w| rings — highlights poles and zeros), "Cyclic", "Rainbow", "CoolTones", "WarmTones", "Greyscale", "Temperature" ColorFunctionScaling True (default): scale each of the eight arguments to \[0,1\] across the sampled domain before calling a custom ColorFunction RegionFunction      f\[x,y\] mask; excluded cells are not drawn PlotLegends         Automatic / True: attach a vertical phase color scale bar (thermal ramp, -π at bottom, π at top) Standard Graphics options (Axes, AspectRatio→1, Frame, PlotRange, AxesLabel, GridLines, ImageSize, Background, PlotLabel, …) pass through to the Graphics\[...\] result. Examples: ComplexPlot\[z^2, {z, -2-2I, 2+2I}\] ComplexPlot\[Sin\[z\], {z, -Pi-Pi\*I, Pi+Pi\*I}\] ComplexPlot\[1/(z^2+1), {z, -2-2I, 2+2I}, PlotPoints-\>80\] ComplexPlot\[(z^2+1)/(z^2-1), {z, -2-2I, 2+2I}, PlotLegends-\>Automatic\]

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

**See also:** [HoldAll](../../expression-information/HoldAll/), [Plot](../../graphics/Plot/), [Arg](../../arithmetic/Arg/), [AspectRatio](../../other-advanced/AspectRatio/), [Frame](../../other-advanced/Frame/), [ImageSize](../../other-advanced/ImageSize/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
