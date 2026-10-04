# Plot3D

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Plot3D[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]`**

Samples f over a uniform grid on \[xmin,xmax\] x \[ymin,ymax\], displays the resulting surface in an interactive orbit-camera window, and returns it as a Graphics3D\[...\] object. A list of functions Plot3D\[{f1, f2, ...}, {x,...}, {y,...}\] draws each surface in a distinct palette colour. Shares Plot's option semantics where they apply: PlotPoints (per-axis grid resolution, default 25), MaxRecursion (doubles the whole grid's resolution while a flatness check fails, default 2 -- a global, crack-free analogue of Plot's adaptive bisection), Mesh (overlay the grid wireframe; default True, unlike Plot's None), PlotStyle, ColorFunction (a function of scaled-x and z, or "Rainbow"), ColorFunctionScaling (default True), RegionFunction (f\[x,y,z\], or Plot's f\[x,y\]/f\[x\] forms), PlotRange (an explicit {zmin,zmax} z-band), Axes, PlotLabel, Background, ImageSize, Lighting (Automatic (default, Lambertian shading) or None to disable).

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Plot3D[{Sin[x + y], Cos[x - y]}, {x, -2, 2}, {y, -2, 2}]
Out[1]= -Graphics-
```

### Options (5)

```mathematica
In[2]:= Plot3D[x^2 - y^2, {x, -2, 2}, {y, -2, 2}, ColorFunction -> "Rainbow", Mesh -> None]
Out[2]= -Graphics-

In[3]:= Plot3D[x + y, {x, -2, 2}, {y, -2, 2}, RegionFunction -> Function[{x, y, z}, x^2 + y^2 < 4]]
Out[3]= -Graphics-

In[4]:= Plot3D[x + y, {x,-2,2}, {y,-2,2},RegionFunction -> Function[{x,y,z}, x^2+y^2 <4],ExclusionStyle -> RGBColor[1, 0.3, 0]]
Out[4]= -Graphics-

In[5]:= Plot3D[{x^2, x^2 + 1}, {x,-2,2}, {y,-2,2},PlotStyle -> {Blue, Red}]
Out[5]= -Graphics-

In[6]:= Plot3D[Sin[x] Cos[2 y], {x, -5, 5}, {y, -5, 5}, BoxRatios -> {1, 1, 1}]
Out[6]= -Graphics-
```

### Applications (5)

```mathematica
In[7]:= Plot3D[Sin[x] Cos[y], {x, -3, 3}, {y, -3, 3}]
Out[7]= -Graphics-

In[8]:= Plot3D[x^2 - y^2, {x, -2, 2}, {y, -2, 2}, ColorFunction -> "Rainbow", Mesh -> None]
Out[8]= -Graphics-

In[9]:= Plot3D[x + y, {x, -2, 2}, {y, -2, 2}, RegionFunction -> Function[{x, y, z}, x^2 + y^2 < 4]]
Out[9]= -Graphics-

In[10]:= Head[Plot3D[x^2 + y^2, {x, -1, 1}, {y, -1, 1}]]
Out[10]= Graphics3D

In[11]:= Options[Plot3D, {PlotPoints, BoxRatios, Mesh}]
Out[11]= {PlotPoints -> 25, BoxRatios -> {1.0, 1.0, 0.4}, Mesh -> Automatic}
```

## Algorithm

plot3d.c — Plot3D[f, {x,xmin,xmax}, {y,ymin,ymax}, opts...].

Mirrors plot.c's shape as closely as the dimensionality allows: HoldAll (the iterator vars have no value yet), a single split_options3() pass that separates the sampler's own options from a passthrough list copied onto the resulting Graphics3D[...] result -- exactly split_options's role in plot.c, just with the smaller set of options that have a 3D meaning. Anything not recognised here falls into the generic passthrough branch and is inertly ignored by the renderer.

The one place 3D genuinely cannot reuse 2D's sampler (sampling.c) is the adaptive refinement itself: that sampler bisects an *ordered* 1D interval, which has no 2D analogue without inventing a per-cell quadtree -- and a quadtree creates T-junction cracks where differently-refined cells meet. Instead MaxRecursion here doubles the *whole* grid's resolution when a cheap flatness spot-check fails, capped at a few levels and a hard point-count ceiling. This stays crack-free (every level is a uniform grid) and gives MaxRecursion real meaning, at the cost of refining more than strictly necessary -- an acceptable trade for "simple and clear" over a true adaptive mesh.

## Implementation notes

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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Plot](../../graphics/Plot/), [HoldAll](../../expression-information/HoldAll/), [ImageSize](../../other-advanced/ImageSize/), [Lighting](../../other-advanced/Lighting/)

- Source: [`src/graphics/plot3d.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/plot3d.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_manipulate.c`](https://github.com/stblake/mathilda/blob/main/tests/test_manipulate.c)
- Tests: [`tests/test_plot3d.c`](https://github.com/stblake/mathilda/blob/main/tests/test_plot3d.c)

## Notes & additional examples

### Notes

`Plot3D` is `HoldAll` and samples a **uniform `PlotPoints x PlotPoints` grid**,
turning each cell into one quad `Polygon`. Refinement is not a quadtree: because
bisecting an ordered interval has no crack-free 2-D analogue, `MaxRecursion`
instead doubles the whole grid resolution (capped at 200/axis) while the surface
is not yet flat, so the mesh always stays crack-free. A cell with any corner where
the function is non-finite is dropped, leaving a hole; a cell straddling a
`RegionFunction` boundary is Sutherland–Hodgman-clipped into a smooth partial
polygon.

The default display box is `BoxRatios -> {1, 1, 0.4}`, so a shallow z-range still
reads as a surface; a raw `Graphics3D[...]` renders at true scale. A single
surface with no `ColorFunction` or `PlotStyle` gets an implicit `"Viridis"` height
gradient. The window is an orbit camera (drag to rotate, scroll to zoom); there is
no re-sampling on rotate because the sampled `(x, y)` domain never changes.
