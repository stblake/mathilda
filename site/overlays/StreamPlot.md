### Worked examples

```mathematica
(* evenly-spaced circular streamlines of a rotation field *)
In[1]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}]
```

```mathematica
(* a nonlinear field *)
In[1]:= StreamPlot[{1 - y^2, x}, {x, -3, 3}, {y, -2, 2}]
```

```mathematica
(* denser seeding draws more, closer lines *)
In[1]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, StreamPoints -> 30]
```

```mathematica
(* StreamAnimate -> True makes the flow particles move in the interactive window *)
In[1]:= StreamPlot[{-y, x}, {x, -2, 2}, {y, -2, 2}, StreamAnimate -> True]
```

```mathematica
(* HoldAll and Protected *)
In[1]:= Attributes[StreamPlot]
```

### Notes

`StreamPlot` is `HoldAll`. It traces **evenly-spaced** streamlines by the
Jobard–Lefebvre scheme: each line is grown in both directions from a seed by RK4
integration of the *normalised* field (a fixed arc-length step, so point spacing and
rendered curvature stay uniform regardless of local speed), and stops when it nears
an existing line, leaves the domain, reaches a critical point, or closes on itself.
A uniform spatial hash whose cell size equals the target line separation keeps the
lines evenly spaced with `O(1)` proximity tests, and seeds are placed outward from
the domain centre.

Each line renders in Mathematica's dashed-arrow style — short curve-following `Line`
dashes capped by filled `Polygon` arrowheads, so the flow direction reads
everywhere. `StreamPoints` sets the density (default 25); `StreamScale` caps a line's
arc length, while the default lets every line run to its natural end.
`StreamAnimate -> True` instead emits one `AnimatedStreamline[...]` per line whose
particle dots flow downstream in the interactive renderer.
