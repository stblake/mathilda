# Plot3D

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Plot3D[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]`**

Samples f over a uniform grid on \[xmin,xmax\] x \[ymin,ymax\], displays the resulting surface in an interactive orbit-camera window, and returns it as a Graphics3D\[...\] object. A list of functions Plot3D\[{f1, f2, ...}, {x,...}, {y,...}\] draws each surface in a distinct palette colour. Shares Plot's option semantics where they apply: PlotPoints (per-axis grid resolution, default 25), MaxRecursion (doubles the whole grid's resolution while a flatness check fails, default 2 -- a global, crack-free analogue of Plot's adaptive bisection), Mesh (overlay the grid wireframe; default True, unlike Plot's None), PlotStyle, ColorFunction (a function of scaled-x and z, or "Rainbow"), ColorFunctionScaling (default True), RegionFunction (f\[x,y,z\], or Plot's f\[x,y\]/f\[x\] forms), PlotRange (an explicit {zmin,zmax} z-band), Axes, PlotLabel, Background, ImageSize, Lighting (Automatic (default, Lambertian shading) or None to disable).

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Plot3D[Sin[x] Cos[y], {x, -3, 3}, {y, -3, 3}]
Out[1]= -Graphics-

In[2]:= Plot3D[{Sin[x + y], Cos[x - y]}, {x, -2, 2}, {y, -2, 2}]
Out[2]= -Graphics-
```

### Options (4)

```mathematica
In[3]:= Plot3D[x^2 - y^2, {x, -2, 2}, {y, -2, 2}, ColorFunction -> "Rainbow", Mesh -> None]
Out[3]= -Graphics-

In[4]:= Plot3D[x + y, {x, -2, 2}, {y, -2, 2}, RegionFunction -> Function[{x, y, z}, x^2 + y^2 < 4]]
Out[4]= -Graphics-

In[5]:= Plot3D[x + y, {x,-2,2}, {y,-2,2},RegionFunction -> Function[{x,y,z}, x^2+y^2 <4],ExclusionStyle -> RGBColor[1, 0.3, 0]]
Out[5]= -Graphics-

In[6]:= Plot3D[{x^2, x^2 + 1}, {x,-2,2}, {y,-2,2},PlotStyle -> {Blue, Red}]
Out[6]= -Graphics-
```

## Algorithm

plot3d.c — Plot3D[f, {x,xmin,xmax}, {y,ymin,ymax}, opts...].

Mirrors plot.c's shape as closely as the dimensionality allows: HoldAll (the iterator vars have no value yet), a single split_options3() pass that separates the sampler's own options from a passthrough list copied onto the resulting Graphics3D[...] result -- exactly split_options's role in plot.c, just with the smaller set of options that have a 3D meaning. Anything not recognised here falls into the generic passthrough branch and is inertly ignored by the renderer.

The one place 3D genuinely cannot reuse 2D's sampler (sampling.c) is the adaptive refinement itself: that sampler bisects an *ordered* 1D interval, which has no 2D analogue without inventing a per-cell quadtree -- and a quadtree creates T-junction cracks where differently-refined cells meet. Instead MaxRecursion here doubles the *whole* grid's resolution when a cheap flatness spot-check fails, capped at a few levels and a hard point-count ceiling. This stays crack-free (every level is a uniform grid) and gives MaxRecursion real meaning, at the cost of refining more than strictly necessary -- an acceptable trade for "simple and clear" over a true adaptive mesh.

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Plot](../../graphics/Plot/), [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_plot3d.c`](https://github.com/stblake/mathilda/blob/main/tests/test_plot3d.c)
