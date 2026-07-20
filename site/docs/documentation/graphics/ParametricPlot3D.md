# ParametricPlot3D

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
ParametricPlot3D[{fx, fy, fz}, {t, tmin, tmax}, opts...]
    Adaptively samples the parametric 3D space curve (fx(t), fy(t), fz(t)) over [tmin, tmax] and returns a Graphics3D[...] object (auto-displayed in an orbit-camera window). The body may be any expression that evaluates to a 3-element {x,y,z} list. Multiple curves: ParametricPlot3D[{{fx1,fy1,fz1}, ...}, {t,...}].
ParametricPlot3D[{fx, fy, fz}, {t, tmin, tmax}, {u, umin, umax}, opts...]
    Two-iterator form: samples a PlotPoints x PlotPoints grid of (t,u) pairs, maps each to {x,y,z}, and emits Polygon[] quads — a parametric 3D surface patch. Options: PlotPoints (initial sample count/grid size, default 25), MaxRecursion (adaptive refinement depth for curves, default 6), MaxPlotPoints (overall point cap, default Infinity), Mesh (All/True: overlays sample dots for curves or grid lines for surfaces; default None), PlotLegends (Automatic/"Expressions"/{labels...}), ColorFunction ("Rainbow" or f[x,y,z] receiving scaled spatial coords, or f[x,z] / f[z] for height-based coloring), ColorFunctionScaling (default True), RegionFunction (f[x,y,z] mask; falls back to f[x,y] forms), PlotStyle, Axes, PlotRange, AxesLabel, PlotLabel, Background, ImageSize (all passed through to Graphics3D). Lighting -> None disables shading (flat colors); default is Automatic (Lambertian shading, same as Plot3D).
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
