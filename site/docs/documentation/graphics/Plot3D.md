# Plot3D

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Plot3D[f, {x, xmin, xmax}, {y, ymin, ymax}, opts...]
    Samples f over a uniform grid on [xmin,xmax] x [ymin,ymax], displays the resulting surface in an interactive orbit-camera window, and returns it as a Graphics3D[...] object. A list of functions Plot3D[{f1, f2, ...}, {x,...}, {y,...}] draws each surface in a distinct palette colour. Shares Plot's option semantics where they apply: PlotPoints (per-axis grid resolution, default 25), MaxRecursion (doubles the whole grid's resolution while a flatness check fails, default 2 -- a global, crack-free analogue of Plot's adaptive bisection), Mesh (overlay the grid wireframe; default True, unlike Plot's None), PlotStyle, ColorFunction (a function of scaled-x and z, or "Rainbow"), ColorFunctionScaling (default True), RegionFunction (f[x,y,z], or Plot's f[x,y]/f[x] forms), PlotRange (an explicit {zmin,zmax} z-band), Axes, PlotLabel, Background, ImageSize, Lighting (Automatic (default, Lambertian shading) or None to disable).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Plot3D[Sin[x] Cos[y], {x, -3, 3}, {y, -3, 3}]
Out[1]= -Graphics3D-

In[2]:= Plot3D[x^2 - y^2, {x, -2, 2}, {y, -2, 2}, ColorFunction -> "Rainbow", Mesh -> None]
Out[2]= -Graphics3D-

In[3]:= Plot3D[x + y, {x, -2, 2}, {y, -2, 2}, RegionFunction -> Function[{x, y, z}, x^2 + y^2 < 4]]
Out[3]= -Graphics3D-

In[4]:= Plot3D[x + y, {x,-2,2}, {y,-2,2},RegionFunction -> Function[{x,y,z}, x^2+y^2 <4],ExclusionStyle -> RGBColor[1, 0.3, 0]]
Out[4]= -Graphics3D-

In[5]:= Plot3D[{Sin[x + y], Cos[x - y]}, {x, -2, 2}, {y, -2, 2}]
Out[5]= -Graphics3D-

In[6]:= Plot3D[{x^2, x^2 + 1}, {x,-2,2}, {y,-2,2},PlotStyle -> {Blue, Red}]
Out[6]= -Graphics3D-
```

## Implementation notes

- `HoldAll`, `Protected`.
- Declines to evaluate if either iterator spec isn't `{var, min, max}` with

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/graphics/graphics_init.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/graphics_init.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
