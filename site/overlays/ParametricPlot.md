### Worked examples

```mathematica
(* a parametric curve: the unit circle, drawn at AspectRatio 1 *)
In[1]:= ParametricPlot[{Cos[t], Sin[t]}, {t, 0, 2 Pi}]
```

```mathematica
(* a Lissajous figure *)
In[1]:= ParametricPlot[{Sin[2 t], Sin[3 t]}, {t, 0, 2 Pi}]
```

```mathematica
(* a computed body (not a literal pair) works too *)
In[1]:= ParametricPlot[2 {Cos[t], Sin[t]}, {t, 0, 2 Pi}]
```

```mathematica
(* two iterators give a filled region: an annulus r in [1, 2] *)
In[1]:= ParametricPlot[{r Cos[t], r Sin[t]}, {t, 0, 2 Pi}, {r, 1, 2}]
```

```mathematica
(* HoldAll and Protected *)
In[1]:= Attributes[ParametricPlot]
```

### Notes

The one-iterator form draws a curve `{fx(t), fy(t)}`; the two-iterator form fills
a region `{fx(t, r), fy(t, r)}` with `Polygon` quads over a uniform grid. The body
is any expression that evaluates to a 2-element list of finite reals — a literal
`{fx, fy}` takes a compiled fast path, a computed form such as `2 {Cos[t],
Sin[t]}` takes the interpreter path.

The curve form runs its **own** adaptive sampler (parallel to `Plot`'s), refining
by Euclidean deviation from the chord normalised by the bounding-box diagonal,
since neither axis is privileged for a parametric curve; invalid samples break the
curve into separate `Line[...]` runs. The region form is a uniform grid with no
refinement (`Mesh -> All` overlays its grid lines). Default `AspectRatio -> 1`.
