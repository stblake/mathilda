# PolarPlot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PolarPlot[r, {theta, tmin, tmax}, opts...]`**

Plots the polar curve r(theta) by converting to Cartesian coordinates {r\*Cos\[theta\], r\*Sin\[theta\]} and sampling adaptively over \[tmin, tmax\]. Returns a Graphics\[...\] object (auto-displayed).

**`PolarPlot[{r1, r2, ...}, {theta, tmin, tmax}, opts...]`**

Multiple polar curves in distinct palette colours. Negative r values are plotted in the opposite direction (standard polar convention). Default AspectRatio -\> 1 (equal axes). Options (same as ParametricPlot): PlotPoints          - initial sample count per curve (default 75) MaxRecursion        - adaptive refinement depth (default 6) MaxPlotPoints       - total point cap (default Infinity) Mesh                - All/True: overlay evaluation dots; None (default) ColorFunction       - f\[t\] or "Rainbow" (sweeps scaled theta) ColorFunctionScaling - True (default): normalise theta to \[0,1\] RegionFunction      - f\[x,y\] mask PlotStyle           - color/style directive(s) PlotLegends         - Automatic / "Expressions" / label list PolarAxes           - option keyword (accepted; polar grid overlay not yet rendered) Standard Graphics options (AspectRatio, Axes, PlotRange, AxesLabel, Frame, GridLines, PlotLabel, Background, ImageSize, Prolog, Epilog) pass through to the Graphics\[...\] result.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= PolarPlot[1, {t, 0, 2Pi}]
Out[1]= -Graphics-

In[2]:= PolarPlot[Sin[2t], {t, 0, 2Pi}]
Out[2]= -Graphics-
```

### Options (3)

```mathematica
In[3]:= PolarPlot[{1, 2}, {t, 0, 2Pi}, PlotStyle -> {Blue, Red}]
Out[3]= -Graphics-

In[4]:= PolarPlot[t, {t, 0, 4Pi}, ColorFunction -> "Rainbow"]
Out[4]= -Graphics-

In[5]:= PolarPlot[Sin[2t], {t, 0, 2Pi}, Mesh -> All, PlotLabel -> "Rose"]
Out[5]= -Graphics-
```

### Applications (5)

```mathematica
In[6]:= PolarPlot[1, {t, 0, 2 Pi}]
Out[6]= -Graphics-

In[7]:= PolarPlot[Sin[2 t], {t, 0, 2 Pi}]
Out[7]= -Graphics-

In[8]:= PolarPlot[t, {t, 0, 4 Pi}, ColorFunction -> "Rainbow"]
Out[8]= -Graphics-

In[9]:= PolarPlot[{1, 2}, {t, 0, 2 Pi}, PlotStyle -> {Blue, Red}]
Out[9]= -Graphics-

In[10]:= Head[PolarPlot[1 + Cos[t], {t, 0, 2 Pi}]]
Out[10]= Graphics
```

## Algorithm

polarplot.c — PolarPlot[r, {theta, tmin, tmax}, opts...]

HoldAll: r and the iterator are unevaluated on entry, exactly like ParametricPlot. The implementation converts the polar body r(theta) into the Cartesian pair {r*Cos[theta], r*Sin[theta]} and delegates to builtin_parametricplot, so all adaptive sampling, option handling, multi-curve paletting, ColorFunction, Mesh, PlotLegends, etc. come for free without duplicating any logic.

## Implementation notes

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

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [ParametricPlot](../../graphics/ParametricPlot/)

- Source: [`src/graphics/polarplot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/polarplot.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)

## Notes & additional examples

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
