---
source: src/graphics/plot.c
---
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
