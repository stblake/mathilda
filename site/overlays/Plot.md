### Worked examples

```mathematica
In[1]:= Plot[Sin[x], {x, 0, 2 Pi}]
```

```mathematica
(* three curves in the ColorData[97] palette *)
In[1]:= Plot[{Sin[x], Cos[x], Sin[2 x]}, {x, 0, 2 Pi}]
```

```mathematica
(* an asymptote: the y-range is clipped to a robust median/MAD band and the curve breaks at each pole *)
In[1]:= Plot[Tan[x], {x, -2 Pi, 2 Pi}]
```

```mathematica
(* the non-real part x < 0 is dropped, so one contiguous Line run survives *)
In[1]:= Length[Cases[Plot[Sqrt[x], {x, -1, 1}], _Line, Infinity]]
```

```mathematica
(* the adaptive-sampler defaults *)
In[1]:= Options[Plot, {PlotPoints, MaxRecursion}]
```

### Notes

`Plot` is `HoldAll`: the body and the `{x, xmin, xmax}` iterator are held until
sampling binds `x`. The curve is **adaptively** sampled — the shared sampler in
`src/graphics/sampling.c` refines each interval by vertical deviation from the
chord (three interior probes, to defeat periodic aliasing), down to
`MaxRecursion` levels. A singularity, a non-real stretch, or a
`RegionFunction` rejection breaks the polyline into separate `Line[...]`
segments rather than bridging the gap, which is why `Sqrt[x]` over a range that
dips below zero yields a single run.

The returned value is an inert `Graphics[{Line[...], ...}, opts...]` object;
the REPL front end renders any top-level `Graphics` automatically, so a
`Graphics` result prints as `-Graphics-`. A hidden `$PlotResample[...]`
metadata node lets the interactive window re-run the sampler over the visible
band when you zoom, so a magnified curve keeps full resolution instead of
exposing the home grid.
