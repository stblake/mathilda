---
references:
  - "E. Wegert, *Visual Complex Functions: An Introduction with Phase Portraits*, Birkhäuser (2012) — domain colouring / phase portraits."
source: src/graphics/complexplot.c
---
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
