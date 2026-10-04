### Worked examples

```mathematica
In[1]:= Plot3D[Sin[x] Cos[y], {x, -3, 3}, {y, -3, 3}]
```

```mathematica
(* z-height maps to hue under a named ramp *)
In[1]:= Plot3D[x^2 - y^2, {x, -2, 2}, {y, -2, 2}, ColorFunction -> "Rainbow", Mesh -> None]
```

```mathematica
(* a RegionFunction leaves a smooth clipped boundary, not a staircase *)
In[1]:= Plot3D[x + y, {x, -2, 2}, {y, -2, 2}, RegionFunction -> Function[{x, y, z}, x^2 + y^2 < 4]]
```

```mathematica
(* the head is an inert Graphics3D[] object *)
In[1]:= Head[Plot3D[x^2 + y^2, {x, -1, 1}, {y, -1, 1}]]
```

```mathematica
(* the grid, box-ratio and mesh defaults *)
In[1]:= Options[Plot3D, {PlotPoints, BoxRatios, Mesh}]
```

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
