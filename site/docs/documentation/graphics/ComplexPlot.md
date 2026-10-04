# ComplexPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ComplexPlot[f, {z, zmin, zmax}, opts...]`**

Domain-colouring plot of the complex function f over the rectangular region in the complex plane with corners zmin and zmax.  z is bound to Complex\[x, y\] at each grid point; f must return a complex or real number.  Each cell's hue encodes Arg(f(z)) on the cyclic "Cyclic" ramp, and |f(z)| sets an HSL lightness: zeros fade to black, poles to white.  ComplexPlot is HoldAll: f and the iterator spec are held unevaluated until z is given a numeric complex value. Options: PlotPoints          grid resolution per axis (default 400) ColorFunction       a named ramp string, or a function of the eight Mathematica arguments Re\[z\], Im\[z\], Abs\[z\], Arg\[z\], Re\[f\], Im\[f\], Abs\[f\], Arg\[f\] (so #8 is the phase of the value).  Named ramps: "PhaseRings" (hue=phase, brightness=log|w| rings — highlights poles and zeros), "Cyclic", "Rainbow", "CoolTones", "WarmTones", "Greyscale", "Temperature" ColorFunctionScaling True (default): scale each of the eight arguments to \[0,1\] across the sampled domain before calling a custom ColorFunction RegionFunction      f\[x,y\] mask; excluded cells are not drawn PlotLegends         Automatic / True: attach a vertical phase color scale bar (thermal ramp, -π at bottom, π at top) Standard Graphics options (Axes, AspectRatio→1, Frame, PlotRange, AxesLabel, GridLines, ImageSize, Background, PlotLabel, …) pass through to the Graphics\[...\] result. Examples: ComplexPlot\[z^2, {z, -2-2I, 2+2I}\] ComplexPlot\[Sin\[z\], {z, -Pi-Pi\*I, Pi+Pi\*I}\] ComplexPlot\[1/(z^2+1), {z, -2-2I, 2+2I}, PlotPoints-\>80\] ComplexPlot\[(z^2+1)/(z^2-1), {z, -2-2I, 2+2I}, PlotLegends-\>Automatic\]

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Options[ComplexPlot, {PlotPoints, ColorFunction}]
Out[1]= {PlotPoints -> 400, ColorFunction -> "Cyclic"}
```

### Options (5)

```mathematica
In[2]:= ComplexPlot[z^2, {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 40]
Out[2]= -Graphics-

In[3]:= ComplexPlot[(z^2 + 1)/(z^2 - 1), {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 40, PlotLegends -> Automatic]
Out[3]= -Graphics-

In[4]:= ComplexPlot[(z^3 - 3)/z, {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 40, ColorFunction -> "PhaseRings"]
Out[4]= -Graphics-

In[5]:= ComplexPlot[Sin[z], {z, -Pi - Pi I, Pi + Pi I}, PlotPoints -> 30, ColorFunction -> (Hue[#8 + 0.5] &)]
Out[5]= ComplexPlot[Sin[z], {z, -Pi - Pi I, Pi + Pi I}, PlotPoints -> 30, ColorFunction -> (Hue[#8 + 0.5] &)]

In[6]:= Head[ComplexPlot[z, {z, -1 - I, 1 + I}, PlotPoints -> 10]]
Out[6]= Graphics
```

### Applications (5)

```mathematica
In[7]:= ComplexPlot[z^2, {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 40]
Out[7]= -Graphics-

In[8]:= ComplexPlot[1/(z^2 + 1), {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 40]
Out[8]= -Graphics-

In[9]:= ComplexPlot[(z^2 + 1)/(z^2 - 1), {z, -2 - 2 I, 2 + 2 I}, PlotPoints -> 40, ColorFunction -> "PhaseRings"]
Out[9]= -Graphics-

In[10]:= ComplexPlot[Sin[z], {z, -Pi - Pi I, Pi + Pi I}, PlotPoints -> 30, ColorFunction -> (Hue[#8 + 0.5] &)]
Out[10]= ComplexPlot[Sin[z], {z, -Pi - Pi I, Pi + Pi I}, PlotPoints -> 30, ColorFunction -> (Hue[#8 + 0.5] &)]

In[11]:= Options[ComplexPlot, {PlotPoints, ColorFunction}]
Out[11]= {PlotPoints -> 400, ColorFunction -> "Cyclic"}
```

## Algorithm

complexplot.c — ComplexPlot and ComplexPlot3D.

Both functions share the same domain-parsing logic and option set. The only structural difference is what they build from the evaluated grid: ComplexPlot emits Rectangle primitives into Graphics[], and ComplexPlot3D emits Polygon quads (height = |w|, colour = arg(w)) into Graphics3D[].

Coloring convention: the default maps the phase arg(w) onto the cyclic "Cyclic" ramp (t = (atan2(im, re) + π) / (2π) ∈ [0, 1]) and folds the modulus in as an HSL lightness — zeros fade to black, poles to white — so ComplexPlot[f] and ComplexPlot[f, ColorFunction -> "Cyclic"] are identical. A custom ColorFunction receives the eight Mathematica arguments

```text
  Re[z], Im[z], Abs[z], Arg[z], Re[f], Im[f], Abs[f], Arg[f]
```

(so #8 is the phase of the value); with ColorFunctionScaling→True (default) each is scaled to [0,1] across the sampled domain.

Both are HoldAll: the body and the iterator spec are held unevaluated until z is bound to Complex[x, y] at each grid point.

Domain spec:

```text
  {z, zmin, zmax}  — z is the complex iterator variable; zmin and zmax
  are evaluated and their Re/Im parts define the rectangular plotting
  domain: xmin=Re(zmin), xmax=Re(zmax), ymin=Im(zmin), ymax=Im(zmax).
  Both endpoints may be real (imaginary part = 0). 
```

## Implementation notes

**Algorithm.** `builtin_complexplot` is `HoldAll` and shares its file (and its
domain-parsing, option-splitting and colour machinery) with
[`ComplexPlot3D`](ComplexPlot3D.md). `parse_complex_iterator`/`bound_to_complex`
read `{z, zmin, zmax}` as the rectangle `xmin = Re(zmin)`, `xmax = Re(zmax)`,
`ymin = Im(zmin)`, `ymax = Im(zmax)` (endpoints may be pure real; a degenerate
edge is rejected). At each grid point `z` is bound to `Complex[x, y]` and the body
evaluated by `cp_eval`, compiled once with the genuinely complex `autocompile_new_z`
(interpreter fallback at compiled poles). It is a **domain-colouring** plot: the
`(N+1) x (N+1)` `CGrid` is drawn cell by cell as a `Rectangle` averaging the four
corner values. The default colouring (`cp_ramp_color`) maps the phase
`t = (atan2(im, re) + Pi)/(2 Pi)` onto the cyclic `"Cyclic"` ramp (no seam at
`±Pi`) and folds the modulus in as an **HSL lightness** `L = |w|/(1 + |w|)` (= ½
at `|w| = 1`): below ½ the hue fades toward black (zeros → black), above it toward
white (poles → white) — so `ComplexPlot[f]` equals `ComplexPlot[f,
ColorFunction -> "Cyclic"]`. A custom `ColorFunction` receives the **eight**
arguments `Re[z], Im[z], Abs[z], Arg[z], Re[f], Im[f], Abs[f], Arg[f]` (so `#8` is
the value's phase) via `cf_eight`, each scaled to `[0,1]` over the grid when
`ColorFunctionScaling -> True`; `"PhaseRings"` is special-cased to one brightness
ring per e-fold of `|w|`. Pole chatter (`Power::infy`) is muted with
`arith_warnings_mute_push/pop` during sampling and such cells are dropped.
`PlotLegends` emits `emit_phase_color_bar` = `$StreamColorBar[-Pi, Pi, cfn]`. The
result is an inert `Graphics[prims, opts...]`.

**Data structures.** `CPlotOpts`; `CGrid {double re, im; bool valid}` of
`(N+1)^2` cells (`build_cgrid`, with `autocompile_new_z`); `CFRange` holding the
per-argument min/max used to scale the eight colour inputs; `Expr** prims` (a
colour directive and a `Rectangle` per cell).

**Complexity / limits.** `O((N+1)^2)` evaluations with `PlotPoints` default
**400** (high, because the cell pitch is the on-screen resolution) — ~160k
samples. Defaults `Frame -> True`, `Axes -> False`, `AspectRatio -> 1`.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [Plot](../../graphics/Plot/), [Arg](../../arithmetic/Arg/), [AspectRatio](../../other-advanced/AspectRatio/), [Frame](../../other-advanced/Frame/), [ImageSize](../../other-advanced/ImageSize/)

- E. Wegert, *Visual Complex Functions: An Introduction with Phase Portraits*, Birkhäuser (2012) — domain colouring / phase portraits.
- Source: [`src/graphics/complexplot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/complexplot.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)

## Notes & additional examples

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
