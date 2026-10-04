### Worked examples

```mathematica
(* domain colouring: hue is the phase, lightness folds in the modulus *)
In[1]:= ComplexPlot[z^2, {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 40]
```

```mathematica
(* poles at +-I fade toward white, the zero at the origin toward black *)
In[1]:= ComplexPlot[1/(z^2 + 1), {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 40]
```

```mathematica
(* "PhaseRings" adds one brightness ring per e-fold of |f| *)
In[1]:= ComplexPlot[(z^2 + 1)/(z^2 - 1), {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 40, ColorFunction -> "PhaseRings"]
```

```mathematica
(* a custom map: #8 is the value's phase Arg[f] *)
In[1]:= ComplexPlot[Sin[z], {z, -Pi - Pi I, Pi + Pi I}, PlotPoints -> 30, ColorFunction -> (Hue[#8 + 0.5] &)]
```

```mathematica
(* the cyclic phase ramp and high grid default *)
In[1]:= Options[ComplexPlot, {PlotPoints, ColorFunction}]
```

### Notes

`ComplexPlot` is `HoldAll`. It binds `z` to `Complex[x, y]` over the rectangle with
corners `zmin`, `zmax` and draws a **domain-colouring** plot: the default maps the
phase `Arg[f(z)]` onto the cyclic `"Cyclic"` ramp (no colour seam at `±Pi`) and
folds the modulus in as an HSL lightness `L = |f|/(1 + |f|)`, so zeros fade to black
and poles to white. `ComplexPlot[f]` is therefore identical to `ComplexPlot[f,
ColorFunction -> "Cyclic"]`.

A custom `ColorFunction` receives **eight** arguments —
`Re[z], Im[z], Abs[z], Arg[z], Re[f], Im[f], Abs[f], Arg[f]` — so `#8` is the
value's phase; `"PhaseRings"` is a built-in that highlights poles and zeros as
concentric ring clusters. The `PlotPoints` default is a high 400 because the cell
pitch is the on-screen resolution; points landing on a pole are dropped silently
(the `Power::infy` chatter is muted during sampling, as in `Plot`).
