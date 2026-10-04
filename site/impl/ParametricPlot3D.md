---
source: src/graphics/parametricplot3d.c
---
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
