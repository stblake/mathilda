# ComplexPlot3D

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ComplexPlot3D[f, {z, zmin, zmax}, opts...]`**

Three-dimensional surface plot of a complex function: height = |f(z)|, colour = Arg(f(z)) via the thermal ramp (same default as Plot3D). z is bound to Complex\[x, y\] at each grid point; the result is a Graphics3D\[...\] object rendered in an orbit-camera window. ComplexPlot3D is HoldAll. Options: PlotPoints          grid resolution per axis (default 200) ColorFunction       a named ramp string, or a function of the eight arguments Re\[z\], Im\[z\], Abs\[z\], Arg\[z\], Re\[f\], Im\[f\], Abs\[f\], Arg\[f\] (see ComplexPlot); "PhaseRings" recommended ColorFunctionScaling True (default): scale the eight args to \[0,1\] RegionFunction      f\[x,y\] mask PlotLegends         Automatic / True: attach a vertical phase color scale bar (-π at bottom, π at top) Lighting -\> None    disables Lambertian shading (recommended for accurate phase colours; default Automatic) Standard Graphics3D options pass through to the result. Examples: ComplexPlot3D\[z^2, {z, -2-2I, 2+2I}\] ComplexPlot3D\[Sin\[z\], {z, -2-2I, 2+2I}, Lighting-\>None\] ComplexPlot3D\[1/z, {z, -2-2I, 2+2I}, PlotLegends-\>Automatic, Lighting-\>None\]

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Options[ComplexPlot3D, {PlotPoints, Lighting}]
Out[1]= {PlotPoints -> 200, Lighting -> Automatic}
```

### Options (3)

```mathematica
In[2]:= ComplexPlot3D[z^2, {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 30]
Out[2]= -Graphics-

In[3]:= ComplexPlot3D[1/z, {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 30, Lighting -> None, PlotLegends -> Automatic]
Out[3]= -Graphics-

In[4]:= Head[ComplexPlot3D[Sin[z], {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 20]]
Out[4]= Graphics3D
```

### Applications (4)

```mathematica
In[5]:= ComplexPlot3D[z^2, {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 30]
Out[5]= -Graphics-

In[6]:= ComplexPlot3D[1/z, {z, -1 - I, 1 + I}, PlotPoints -> 30, Lighting -> None, PlotLegends -> Automatic]
Out[6]= -Graphics-

In[7]:= Head[ComplexPlot3D[Sin[z], {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 20]]
Out[7]= Graphics3D

In[8]:= Options[ComplexPlot3D, {PlotPoints, Lighting}]
Out[8]= {PlotPoints -> 200, Lighting -> Automatic}
```

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

**Algorithm.** `builtin_complexplot3d` is `HoldAll` and shares the same file and
substrate as [`ComplexPlot`](ComplexPlot.md): identical `{z, zmin, zmax}` domain
parsing (`parse_complex_iterator`/`bound_to_complex`), the same
`split_cplot_options`/`CPlotOpts` option set, the same `cp_eval` body evaluation
(compiled with `autocompile_new_z`), the same `CGrid`, and the same eight-argument
`ColorFunction`/`ColorFunctionScaling` convention. The difference is that it
builds a **surface**: per vertex the **height is `|f(z)|`** and the default
**colour is `Arg(f(z))` on the thermal ramp** (`thermal_rgb`) — crucially
*without* the modulus-brightness attenuation used by the 2-D head, so the surface
keeps the same visual weight as `Plot3D`. Each cell becomes a
`Polygon[{p00, p10, p11, p01}]` quad of `{x, y, z} = {Re(z), Im(z), |f(z)|}`
vertices. To keep a pole from spiking to infinity, `compute_height_cap` takes the
**95th percentile** of the valid `|f|` values as a cap and clamps heights to it,
so a pole renders as a flat-topped column; `embed_plot_range3` supplies the
`{0, hcap}` z-extent. `Lighting` (and any option not consumed by
`split_cplot_options`) rides the generic pass-through onto the `Graphics3D`;
`Lighting -> None` disables Lambertian shading for accurate phase colours. The
result is an inert `Graphics3D[prims, opts...]`.

**Data structures.** Shared with `ComplexPlot`: `CPlotOpts`, `CGrid`, `CFRange`;
plus the sorted valid-`|f|` array used by `compute_height_cap`; `Expr** prims`
(a colour directive and a `Polygon` quad per cell).

**Complexity / limits.** `O((N+1)^2)` evaluations and quads with `PlotPoints`
default **200** — kept well below `ComplexPlot`'s 400 because this is an `N x N`
polygon mesh, not a raster (400² would be ~320k triangles). Defaults
`Frame -> True`, `Axes -> False`, `AspectRatio -> 1`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Plot3D](../../graphics/Plot3D/), [HoldAll](../../expression-information/HoldAll/), [ComplexPlot](../../graphics/ComplexPlot/), [Lighting](../../other-advanced/Lighting/)

- E. Wegert, *Visual Complex Functions: An Introduction with Phase Portraits*, Birkhäuser (2012) — phase-portrait surfaces.
- Source: [`src/graphics/complexplot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/complexplot.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)

## Notes & additional examples

### Notes

`ComplexPlot3D` is `HoldAll` and shares its domain parsing, option set and body
evaluation with [`ComplexPlot`](ComplexPlot.md); the eight-argument `ColorFunction`
convention is identical. The difference is a **surface**: each vertex has height
`|f(z)|` and default colour `Arg(f(z))` on the thermal ramp — with no
modulus-brightness attenuation, so the surface keeps the same visual weight as
`Plot3D`.

A pole would spike to infinity, so heights are clamped to the 95th percentile of
the sampled `|f|` values, rendering a pole as a flat-topped column. `Lighting ->
None` disables Lambertian shading for accurate phase colours. `PlotPoints` defaults
to 200 (well below `ComplexPlot`'s 400) because the surface is an `N x N` polygon
mesh, not a raster.
