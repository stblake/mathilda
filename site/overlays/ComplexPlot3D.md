### Worked examples

```mathematica
(* height = |f(z)|, colour = Arg[f(z)] on the thermal ramp *)
In[1]:= ComplexPlot3D[z^2, {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 30]
```

```mathematica
(* a simple pole, flat-shaded for true phase colours, with the phase legend *)
In[1]:= ComplexPlot3D[1/z, {z, -1 - I, 1 + I}, PlotPoints -> 30, Lighting -> None, PlotLegends -> Automatic]
```

```mathematica
(* the result is an inert Graphics3D[] object *)
In[1]:= Head[ComplexPlot3D[Sin[z], {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 20]]
```

```mathematica
(* a lower grid default than ComplexPlot: this is a polygon mesh, not a raster *)
In[1]:= Options[ComplexPlot3D, {PlotPoints, Lighting}]
```

### Notes

`ComplexPlot3D` is `HoldAll` and shares its domain parsing, option set and body
evaluation with [`ComplexPlot`](ComplexPlot.md); the eight-argument `ColorFunction`
convention is identical. The difference is a **surface**: each vertex has height
`|f(z)|` and default colour `Arg(f(z))` on the thermal ramp — with no
modulus-brightness attenuation, so the surface keeps the same visual weight as
`Plot3D`.

A pole would spike to infinity, so heights are clamped to the 95th percentile of
the sampled `|f|` values, rendering a pole as a flat-topped column. `Lighting ->
None` disables Lambertian shading for accurate phase colours. `PlotPoints` defaults
to 200 (well below `ComplexPlot`'s 400) because the surface is an `N x N` polygon
mesh, not a raster.
