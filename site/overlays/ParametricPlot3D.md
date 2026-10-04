### Worked examples

```mathematica
(* a space curve: a helix *)
In[1]:= ParametricPlot3D[{Cos[t], Sin[t], t/5}, {t, 0, 4 Pi}]
```

```mathematica
(* two iterators give a surface patch: the unit sphere *)
In[1]:= ParametricPlot3D[{Cos[u] Sin[v], Sin[u] Sin[v], Cos[v]}, {u, 0, 2 Pi}, {v, 0, Pi}]
```

```mathematica
(* a torus, R = 2, r = 1 *)
In[1]:= ParametricPlot3D[{(2 + Cos[v]) Cos[u], (2 + Cos[v]) Sin[u], Sin[v]}, {u, 0, 2 Pi}, {v, 0, 2 Pi}]
```

```mathematica
(* the result is an inert Graphics3D[] object *)
In[1]:= Head[ParametricPlot3D[{t, t^2, t^3}, {t, 0, 1}]]
```

```mathematica
(* HoldAll and Protected *)
In[1]:= Attributes[ParametricPlot3D]
```

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
