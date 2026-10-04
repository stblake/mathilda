---
references:
  - "E. Wegert, *Visual Complex Functions: An Introduction with Phase Portraits*, Birkhäuser (2012) — phase-portrait surfaces."
source: src/graphics/complexplot.c
---
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
