---
source: src/graphics/plot3d.c
---
**Algorithm.** `builtin_plot3d` is `HoldAll`. It parses two iterator specs
`{x, xmin, xmax}`, `{y, ymin, ymax}` (both must be `ITER_KIND_RANGE` with ordered
numeric bounds, else `NULL`), splits options with `split_options3` into a
`Plot3DSampleOpts` bundle and an evaluated pass-through array, and accepts a
single surface or a `List` of surfaces. `build_surface_primitives` samples each
surface on a **uniform `n x n` grid** (`n = PlotPoints`) with the body
auto-compiled once by `autocompile_new(body, {x, y}, 2)` (interpreter fallback
per point). Refinement is *not* a quadtree: an ordered 1-D bisection has no
crack-free 2-D analogue, so while `surface_is_flat` fails — the mid-cell value
against the bilinear interpolant of the four corners, to `FLAT_TOL3D = 0.0025` of
the robust z-span — the **whole grid doubles** (`n *= 2`), up to `MaxRecursion`
levels and a hard cap of 200 points/axis, staying uniform (hence crack-free) at
the cost of global over-refinement. Each fully valid cell becomes one 4-vertex
`Polygon[{p00, p10, p11, p01}]`; a cell with any corner where `f` is non-finite
is skipped entirely, leaving a **hole**; a cell straddling a `RegionFunction`
boundary is clipped by Sutherland–Hodgman (`clip_quad_to_region`, crossings found
by bisection) into a smooth partial polygon rather than a staircase. A single
surface with no `ColorFunction`/`PlotStyle` gets an implicit `"Viridis"` height
gradient; multi-surface uses `palette_color`; an explicit `PlotStyle` colour
renders solid. The result is an inert `Graphics3D[prims, opts...]`.

**Data structures.** `Plot3DSampleOpts`; `Plot3DEvalCtx` (`varx`, `vary`,
`body`, `RegionFunction`, `AutoCompiled*`); `GridPt {x, y, z; valid; fn_ok}`
splitting "function defined" from "region-accepted"; `ClipPt {x, y, z;
is_crossing}` for the boundary-clipped polygons. There is no `$PlotResample`
analogue — the 3-D window orbits a fixed mesh rather than re-sampling.

**Complexity / limits.** `O(n^2)` evaluations per refinement level, `n` doubling
to at most 200/axis, plus one centre evaluation per cell for the flatness test;
default grid `PlotPoints = 25`, `MaxRecursion = 2`, `Mesh -> True`,
`BoxRatios -> {1, 1, 0.4}`, `Axes -> True`. The global refinement means one
locally wiggly patch forces the entire grid finer.
