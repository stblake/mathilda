---
source: src/graphics/polarplot.c
---
**Algorithm.** `builtin_polarplot` is a thin `HoldAll` wrapper with no sampling
logic of its own. `make_polar_pair` rewrites the polar body `r` and iterator
variable `theta` into the Cartesian pair `{r Cos[theta], r Sin[theta]}`; a
`List` body `{r1, r2, ...}` maps each `ri` to its own pair, producing a
`List`-of-`List`s so `ParametricPlot`'s multi-curve detection fires. It then
synthesises a `ParametricPlot[param_body, iter, opts..., PlotPoints -> 75]` call
and delegates to `builtin_parametricplot`, freeing the synthetic call afterward.
All adaptive sampling, `ColorFunction`, `Mesh`, `PlotLegends` and option handling
therefore come from `ParametricPlot`; negative `r` falls out naturally from the
`r Cos`/`r Sin` rewrite (the point reflects through the origin). `PolarPlot`
injects `PlotPoints -> 75` only when the caller has not already passed
`PlotPoints` — more initial seeds than `ParametricPlot`'s 50 because a polar
curve spans a full `2 Pi` and tight petals otherwise render angular. The result
is the `Graphics[prims, opts...]` object `ParametricPlot` returns.

**Data structures.** None beyond the synthesised `ParametricPlot[...]` `Expr`
tree handed to the delegate; all sample buffers belong to `ParametricPlot`.

**Complexity / limits.** Same cost profile as `ParametricPlot`'s curve form with
`PlotPoints = 75`, `MaxRecursion = 6`, default `AspectRatio -> 1`, `PlotStyle`
colour `RGBColor[0.2, 0.4, 0.8]`, all inherited. The `PolarAxes` option is
accepted but its polar grid overlay is not yet rendered (see
[`PolarAxes`](PolarAxes.md)).
