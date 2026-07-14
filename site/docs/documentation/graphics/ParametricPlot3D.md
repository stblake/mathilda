# ParametricPlot3D

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ParametricPlot3D[{fx, fy, fz}, {t, tmin, tmax}, opts...]
    Adaptively samples the parametric 3D space curve (fx(t), fy(t), fz(t)) over [tmin, tmax] and returns a Graphics3D[...] object (auto-displayed in an orbit-camera window). The body may be any expression that evaluates to a 3-element {x,y,z} list. Multiple curves: ParametricPlot3D[{{fx1,fy1,fz1}, ...}, {t,...}].
ParametricPlot3D[{fx, fy, fz}, {t, tmin, tmax}, {u, umin, umax}, opts...]
    Two-iterator form: samples a PlotPoints x PlotPoints grid of (t,u) pairs, maps each to {x,y,z}, and emits Polygon[] quads — a parametric 3D surface patch. Options: PlotPoints (initial sample count/grid size, default 25), MaxRecursion (adaptive refinement depth for curves, default 6), MaxPlotPoints (overall point cap, default Infinity), Mesh (All/True: overlays sample dots for curves or grid lines for surfaces; default None), PlotLegends (Automatic/"Expressions"/{labels...}), ColorFunction ("Rainbow" or f[x,y,z] receiving scaled spatial coords, or f[x,z] / f[z] for height-based coloring), ColorFunctionScaling (default True), RegionFunction (f[x,y,z] mask; falls back to f[x,y] forms), PlotStyle, Axes, PlotRange, AxesLabel, PlotLabel, Background, ImageSize (all passed through to Graphics3D). Lighting -> None disables shading (flat colors); default is Automatic (Lambertian shading, same as Plot3D).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ParametricPlot3D[{Cos[t], Sin[t], t/5}, {t, 0, 4 Pi}]
Out[1]= -Graphics3D-

In[2]:= ParametricPlot3D[{Sin[2 t] Cos[t], Sin[2 t] Sin[t], Cos[t]}, {t, 0, 2 Pi}]
Out[2]= -Graphics3D-

In[3]:= ParametricPlot3D[{{Cos[t], Sin[t], 0}, {Cos[t], 0, Sin[t]}}, {t, 0, 2 Pi}]
Out[3]= -Graphics3D-

In[4]:= ParametricPlot3D[{Cos[t], Sin[t], t/5}, {t, 0, 4 Pi}, Mesh -> All]
Out[4]= -Graphics3D-
```

## Implementation notes

- `HoldAll`, `Protected`.
- Declines to evaluate if bounds aren't numeric or `PlotPoints < 2`.
- Auto-displays exactly like `Graphics3D`/`Plot3D`.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
