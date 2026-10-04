# ParametricPlot3D

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ParametricPlot3D[{fx, fy, fz}, {t, tmin, tmax}, opts...]`**

Adaptively samples the parametric 3D space curve (fx(t), fy(t), fz(t)) over \[tmin, tmax\] and returns a Graphics3D\[...\] object (auto-displayed in an orbit-camera window). The body may be any expression that evaluates to a 3-element {x,y,z} list. Multiple curves: ParametricPlot3D\[{{fx1,fy1,fz1}, ...}, {t,...}\].

**`ParametricPlot3D[{fx, fy, fz}, {t, tmin, tmax}, {u, umin, umax}, opts...]`**

Two-iterator form: samples a PlotPoints x PlotPoints grid of (t,u) pairs, maps each to {x,y,z}, and emits Polygon\[\] quads — a parametric 3D surface patch. Options: PlotPoints (initial sample count/grid size, default 25), MaxRecursion (adaptive refinement depth for curves, default 6), MaxPlotPoints (overall point cap, default Infinity), Mesh (All/True: overlays sample dots for curves or grid lines for surfaces; default None), PlotLegends (Automatic/"Expressions"/{labels...}), ColorFunction ("Rainbow" or f\[x,y,z\] receiving scaled spatial coords, or f\[x,z\] / f\[z\] for height-based coloring), ColorFunctionScaling (default True), RegionFunction (f\[x,y,z\] mask; falls back to f\[x,y\] forms), PlotStyle, Axes, PlotRange, AxesLabel, PlotLabel, Background, ImageSize (all passed through to Graphics3D). Lighting -\> None disables shading (flat colors); default is Automatic (Lambertian shading, same as Plot3D).

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= ParametricPlot3D[{Cos[t], Sin[t], t/5}, {t, 0, 4 Pi}]
Out[1]= -Graphics-

In[2]:= ParametricPlot3D[{Sin[2 t] Cos[t], Sin[2 t] Sin[t], Cos[t]}, {t, 0, 2 Pi}]
Out[2]= -Graphics-

In[3]:= ParametricPlot3D[{{Cos[t], Sin[t], 0}, {Cos[t], 0, Sin[t]}}, {t, 0, 2 Pi}]
Out[3]= -Graphics-

In[4]:= ParametricPlot3D[{Cos[u] Sin[v], Sin[u] Sin[v], Cos[v]}, {u, 0, 2 Pi}, {v, 0, Pi}]
Out[4]= -Graphics-

In[5]:= ParametricPlot3D[{(2 + Cos[v]) Cos[u], (2 + Cos[v]) Sin[u], Sin[v]}, {u, 0, 2 Pi}, {v, 0, 2 Pi}]
Out[5]= -Graphics-

In[6]:= ParametricPlot3D[{{Cos[u] Sin[v], Sin[u] Sin[v], Cos[v]}, {2 Cos[u] Sin[v], 2 Sin[u] Sin[v], 2 Cos[v]}}, {u, 0, 2 Pi}, {v, 0, Pi}]
Out[6]= -Graphics-
```

### Options (4)

```mathematica
In[7]:= ParametricPlot3D[{Cos[t], Sin[t], t/5}, {t, 0, 4 Pi}, ColorFunction -> "Rainbow"]
Out[7]= -Graphics-

In[8]:= ParametricPlot3D[{Cos[t], Sin[t], t/5}, {t, 0, 4 Pi}, Mesh -> All]
Out[8]= -Graphics-

In[9]:= ParametricPlot3D[{u Cos[v], u Sin[v], u}, {u, 0, 2}, {v, 0, 2 Pi}, ColorFunction -> "Rainbow", Mesh -> All]
Out[9]= -Graphics-

In[10]:= ParametricPlot3D[{Cos[u] Sin[v], Sin[u] Sin[v], Cos[v]}, {u, 0, 2 Pi}, {v, 0, Pi}, Lighting -> None]
Out[10]= -Graphics-
```

### Applications (5)

```mathematica
In[11]:= ParametricPlot3D[{Cos[t], Sin[t], t/5}, {t, 0, 4 Pi}]
Out[11]= -Graphics-

In[12]:= ParametricPlot3D[{Cos[u] Sin[v], Sin[u] Sin[v], Cos[v]}, {u, 0, 2 Pi}, {v, 0, Pi}]
Out[12]= -Graphics-

In[13]:= ParametricPlot3D[{(2 + Cos[v]) Cos[u], (2 + Cos[v]) Sin[u], Sin[v]}, {u, 0, 2 Pi}, {v, 0, 2 Pi}]
Out[13]= -Graphics-

In[14]:= Head[ParametricPlot3D[{t, t^2, t^3}, {t, 0, 1}]]
Out[14]= Graphics3D

In[15]:= Attributes[ParametricPlot3D]
Out[15]= {HoldAll, Protected}
```

## Algorithm

parametricplot3d.c — ParametricPlot3D[body, {t, tmin, tmax}, opts...]

```text
                    — ParametricPlot3D[body, {t, tmin, tmax}, {u, umin, umax}, opts...]
```

HoldAll: the body and all iterator specs are unevaluated when received.

One-iterator form:

```text
  body must evaluate to {x, y, z} for each t. Single body or
  {body1, body2,...} (list of bodies) for multi-curve.
  Adaptive 3D sampling with a Euclidean chord-deviation test in (x,y,z) space.
  Three interior probes per interval prevent aliasing against periodic curves.
  Output: Graphics3D[{Line[...], ...}, opts].
```

Two-iterator form:

```text
  body evaluates to {x, y, z} for each (t, u) pair.
  Samples a PlotPoints x PlotPoints grid and builds Polygon[] quads,
  like Plot3D but driven by a parametric mapping instead of f(x,y).
  Output: Graphics3D[{Polygon[...], ...}, opts].
```

In both forms the body can be any expression that evaluates to a 3-element numeric list: a literal {fx, fy, fz}, or a computed form, as long as the result has head List and exactly three finite-real elements. ColorFunction receives scaled (x,y,z) coordinates (spatial, not parameter), identical to Plot3D's convention. "Rainbow" sweeps hue over the z-extent.

## Implementation notes

**Algorithm.** `builtin_parametricplot3d` is `HoldAll` and mirrors
`ParametricPlot` one dimension up: one iterator is a **space curve**, two is a
**surface patch**. The body must evaluate to a 3-element list of finite reals;
`Param3DEvalCtx` holds an `ac[3]` compiled triple (`param3d_ctx_compile`,
all-or-nothing on a literal 3-element `List`) with the interpreter fallback
`eval_body_xyz`. `RegionFunction` is evaluated by `eval_region3d`, which tries
`f[x, y, z]` first and falls back to `f[x, y]` when that is non-Boolean. The
curve form (`build_param3d_curve` → `param3d_sample`) uses the same three-probe
adaptive bisection with a Euclidean deviation test in `(x, y, z)` normalised by
the 3-D bounding-box diagonal, at the looser tolerance `PARAM3D_FLAT_TOL =
0.0025`; it emits `Line[...]` runs broken at invalid samples. The surface form
(`build_param3d_surface`) samples a uniform `n x n` `(t, u)` grid and emits one
`Polygon[{p00, p10, p11, p01}]` per all-valid cell, with no adaptivity.
`ColorFunction` receives **scaled spatial** coordinates `{xs, ys, zs}` (not the
parameters), scaled over the sampled xyz bounding box — the same convention as
`Plot3D`, so `"Rainbow"` sweeps hue over the z-extent. A `List`-of-`List`s body
is the multi-object form, each in a plain `palette_color`. The result is an inert
`Graphics3D[prims, opts...]`.

**Data structures.** `Param3DEvalCtx` (iterator vars, body, `ac[3]`,
`RegionFunction`); growing sample buffers for the curve; the option bundle from
`split_options_param3d`.

**Complexity / limits.** `PlotPoints` default **25** for both forms,
`MaxRecursion = 6` for the curve; the surface grid is uniform (a cell with any
invalid corner is dropped). `AspectRatio` is silently discarded (meaningless for
the orbit camera); default `Axes -> True`, `PlotStyle` colour
`RGBColor[0.4, 0.7, 1.0]`. Unlike the 2-D head, multi-object 3-D plots use plain
palette colours with no per-object style scoping.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Plot3D](../../graphics/Plot3D/)

- Source: [`src/graphics/parametricplot3d.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/parametricplot3d.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)

## Notes & additional examples

### Notes

The one-iterator form draws a space curve `{fx(t), fy(t), fz(t)}` with the same
three-probe adaptive sampler as `ParametricPlot`, measuring deviation in
`(x, y, z)` against the 3-D bounding-box diagonal; the two-iterator form tiles a
surface patch `{fx(t, u), fy(t, u), fz(t, u)}` with `Polygon` quads over a uniform
grid. `RegionFunction` is tried as `f[x, y, z]` first, then `f[x, y]`.

`ColorFunction` receives **scaled spatial** coordinates `{xs, ys, zs}` (not the
parameters), the same convention as `Plot3D`, so `"Rainbow"` sweeps hue over the
z-extent. `AspectRatio` is silently dropped (it has no meaning for the orbit
camera), and the window orbits with drag-to-rotate / scroll-to-zoom.
