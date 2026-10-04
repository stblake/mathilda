### Worked examples

```mathematica
(* the unit circle, r = 1 *)
In[1]:= PolarPlot[1, {t, 0, 2 Pi}]
```

```mathematica
(* a four-petal rose *)
In[1]:= PolarPlot[Sin[2 t], {t, 0, 2 Pi}]
```

```mathematica
(* an Archimedean spiral, rainbow-coloured along the curve *)
In[1]:= PolarPlot[t, {t, 0, 4 Pi}, ColorFunction -> "Rainbow"]
```

```mathematica
(* two concentric circles in explicit colours *)
In[1]:= PolarPlot[{1, 2}, {t, 0, 2 Pi}, PlotStyle -> {Blue, Red}]
```

```mathematica
(* a cardioid rewrites to a parametric curve, returning a Graphics[] object *)
In[1]:= Head[PolarPlot[1 + Cos[t], {t, 0, 2 Pi}]]
```

### Notes

`PolarPlot` is a thin `HoldAll` wrapper: it rewrites `r(theta)` into the Cartesian
body `{r Cos[theta], r Sin[theta]}` (mapping each element of a list body to its
own pair) and delegates to `ParametricPlot` with `PlotPoints -> 75`. Everything —
the adaptive sampler, `ColorFunction`, `Mesh`, `PlotLegends`, multi-curve palette
colours — therefore comes from `ParametricPlot`, and a negative `r` falls out
naturally (the point reflects through the origin).

The `PlotPoints` default of 75 (vs `ParametricPlot`'s 50) gives a polar curve
spanning a full `2 Pi` enough initial seeds that tight petals do not render
angular. The `PolarAxes -> True` option is accepted but its polar grid overlay is
not yet drawn — see [`PolarAxes`](PolarAxes.md).
