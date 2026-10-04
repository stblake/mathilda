# Plot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Plot[f, {x, xmin, xmax}, opts...]`**

Adaptively samples f over \[xmin, xmax\], displays the resulting curve in an interactive window, and returns it as a Graphics\[...\] object. A list of functions Plot\[{f1, f2, ...}, {x, xmin, xmax}\] draws each on the same axes in a distinct palette colour. Options: PlotPoints (initial sample count, default 50), MaxRecursion (adaptive refinement depth, default 6), MaxPlotPoints (overall point cap, default Infinity), Mesh (All overlays the evaluation points as dots; default None), PlotRange, PlotRangePadding, AspectRatio, PlotStyle, Axes, AxesLabel, AxesOrigin, AxesStyle, TicksStyle, LabelStyle, Frame, FrameLabel, FrameStyle, FrameTicks, RotateLabel, GridLines, GridLinesStyle, Prolog, Epilog, PlotLabel, Background, ImageSize, ColorFunction (a function, or named ramp: "Rainbow"/"CoolTones"/"WarmTones"/"Greyscale"/"Temperature"), ColorFunctionScaling (default True), Filling (Axis/Bottom/Top/a number), FillingStyle, PlotLegends (Automatic/"Expressions"/an explicit list), RegionFunction, Exclusions.

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Plot[Sin[x], {x, 0, 2 Pi}]
Out[1]= -Graphics-

In[2]:= Plot[1/x, {x, -2, 2}]
Out[2]= -Graphics-

In[3]:= Plot[Sin[x], {x, a, b}]
Out[3]= Plot[Sin[x], {x, a, b}]

In[4]:= Plot[{Sin[1/x], Cos[1/x]}, {x, -Pi, Pi}]
Out[4]= -Graphics-
```

### Options (4)

```mathematica
In[5]:= Plot[Sin[x] + Sin[7 x], {x, -2, 2}, Mesh -> All]
Out[5]= -Graphics-

In[6]:= Plot[Sin[x], {x, 0, 2 Pi}, GridLines -> Automatic, Epilog -> {Red, Point[{0, 0}]}, AxesOrigin -> {0, 0}]
Out[6]= -Graphics-

In[7]:= Plot[Sin[x], {x, 0, 2 Pi}, ColorFunction -> "Rainbow", Filling -> Axis]
Out[7]= -Graphics-

In[8]:= Plot[{Sin[x], Cos[x]}, {x, 0, 2 Pi}, PlotLegends -> "Expressions"]
Out[8]= -Graphics-
```

### Applications (5)

```mathematica
In[9]:= Plot[Sin[x], {x, 0, 2 Pi}]
Out[9]= -Graphics-

In[10]:= Plot[{Sin[x], Cos[x], Sin[2 x]}, {x, 0, 2 Pi}]
Out[10]= -Graphics-

In[11]:= Plot[Tan[x], {x, -2 Pi, 2 Pi}]
Out[11]= -Graphics-

In[12]:= Length[Cases[Plot[Sqrt[x], {x, -1, 1}], _Line, Infinity]]
Out[12]= 1

In[13]:= Options[Plot, {PlotPoints, MaxRecursion}]
Out[13]= {PlotPoints -> 50, MaxRecursion -> 6}
```

## Algorithm

plot.c — Plot[f, {x, xmin, xmax}, opts...].

HoldAll, like Table/Do: f and the iterator spec must not be pre-evaluated (x has no value yet). Splits trailing options into the sampler's own (PlotPoints/MaxRecursion/MaxPlotPoints/Mesh/RegionFunction/ Exclusions/ColorFunction/Filling/..., consumed here) and everything else (PlotRange/AspectRatio/PlotStyle/Axes/.../ImageSize, copied through onto the resulting Graphics[...] unevaluated -- render.c is the single place that interprets those, whether reached via Plot's auto-display or a later Show[]).

## Implementation notes

**Algorithm.** `builtin_plot` is `HoldAll`, so the body `f` and the iterator
`{x, xmin, xmax}` reach it unevaluated. It parses the iterator with
`iter_spec_parse` (requiring an `ITER_KIND_RANGE`), numericises the bounds
through `numericize_bound` (= `expr_to_real_double(N[...])`, so `2 Pi` or
`Sqrt[2]` are valid limits), and declines (returns `NULL`) unless
`xmin < xmax`. `split_options` partitions the trailing `Rule`s into a
`PlotSampleOpts` bundle consumed here and a pass-through array that is each
`evaluate`d once (so `Red` becomes `RGBColor[...]`) and copied onto the result;
a trailing argument that is not a `Rule`, or a known option with a bad value,
makes the whole call decline. A `List` first argument is the multi-curve form.
`build_plot_primitives` shadows the iterator symbol and samples each curve with
`sample_lines`, which: splits the range at `Exclusions`; compiles the body once
with `autocompile_new` (falling back to the interpreter per-point when it
cannot); samples in scaled (`ScalingFunctions`) world space; and feeds each
sub-range to the shared adaptive sampler `plot_sample_adaptive`
(`src/graphics/sampling.c`). That sampler refines **by vertical deviation**: at
each interval midpoint plus the ¼ and ¾ probes it measures the gap to the chord
(clamped into a robust median/MAD display band) and keeps subdividing while
`dev >= FLAT_TOL * yspan` (`FLAT_TOL = 0.0006`) or the chord exceeds
`MAX_CHORD_FRAC = 0.08`, down to `MaxRecursion` levels. Three probes rather than
one defeat periodic aliasing (a `Sin[22 x]` that looks flat on the coarse grid),
and a probe landing in a gap localises a hidden singularity. Non-finite,
complex, or `RegionFunction`-rejected points mark a break, so each contiguous run
of valid points becomes its **own** `Line[...]` — this is how a pole or a
non-real stretch splits the curve instead of bridging it. `ColorFunction` emits
a coloured two-point `Line` per segment; `Filling` emits `Opacity`+`Polygon`
strips; `Mesh` appends the sample `Point`s. The result is an inert
`Graphics[prims, opts..., $PlotResample[...]]`, the last a held metadata node
that lets the renderer re-sample on zoom.

**Data structures.** `PlotSampleOpts` (the full option bundle); `PlotEvalCtx`
(`var`, `body`, `RegionFunction`, the two `ScalingFunctions`, the compiled
`AutoCompiled*`); `PlotPoint {double x, y; bool valid; bool break_before}` from
`sampling.h`; `Range1D` for exclusion sub-ranges. Colours for multi-curve plots
come from `palette_color` — a baked 10-entry `ColorData[97]` table cycled by
`i % 10`. `$PlotResample[var, {bodies}, {opts[10]}]` carries the 10 sampler
parameters so a magnified view re-runs the sampler over the visible band.

**Complexity / limits.** `plot_points` initial evaluations (default **50**) plus
adaptive refinement up to `2^MaxRecursion` (default **6**) extra intervals, ~3–4
evaluations per subdivision, capped by `MaxPlotPoints` (default unbounded).
Default `AspectRatio` is `1/GoldenRatio`, default `PlotStyle` colour
`RGBColor[0.2, 0.4, 0.8]`, `Axes -> True` injected unless `Frame` is given.
Sampling is one-dimensional: `RegionFunction` only rejects points, it does not
clip the curve.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [Show](../../graphics/Show/), [HoldAll](../../expression-information/HoldAll/)

- Source: [`src/graphics/plot.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/plot.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_graphics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics.c)
- Tests: [`tests/test_graphics_sampling.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graphics_sampling.c)
- Tests: [`tests/test_manipulate.c`](https://github.com/stblake/mathilda/blob/main/tests/test_manipulate.c)

## Notes & additional examples

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
