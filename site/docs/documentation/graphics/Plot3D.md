# Plot3D

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Plot3D[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
    Samples f over a uniform grid on [xmin,xmax] x [ymin,ymax], displays the resulting surface in an interactive orbit-camera window, and returns it as a Graphics3D[...] object. A list of functions Plot3D[{f1, f2, ...}, {x,...}, {y,...}] draws each surface in a distinct palette colour. Shares Plot's option semantics where they apply: PlotPoints (per-axis grid resolution, default 25), MaxRecursion (doubles the whole grid's resolution while a flatness check fails, default 2 -- a global, crack-free analogue of Plot's adaptive bisection), Mesh (overlay the grid wireframe; default True, unlike Plot's None), PlotStyle, ColorFunction (a function of scaled-x and z, or "Rainbow"), ColorFunctionScaling (default True), RegionFunction (f[x,y,z], or Plot's f[x,y]/f[x] forms), PlotRange (an explicit {zmin,zmax} z-band), Axes, PlotLabel, Background, ImageSize, Lighting (Automatic (default, Lambertian shading) or None to disable).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
