---
references:
  - "H. A. Sturges, *The Choice of a Class Interval*, J. Amer. Statist. Assoc. **21** (1926) 65-66 — the default bin count."
source: src/graphics/barchart.c
---
**Algorithm.** `builtin_histogram` lives alongside `BarChart` in `barchart.c` and
shares its `split_chart_options` and `emit_bar` machinery. It is a plain
`Protected` builtin (**not** `HoldAll`): the data is evaluated and coerced to
reals. `histogram_parse_bins` reads the bin specification —
`Histogram[data]` uses **Sturges' rule** `n_bins = ceil(log2 n) + 1` clamped to
`[2, 50]`; `Histogram[data, k]` gives `k` bins; `Histogram[data, {step}]` gives
`ceil((dmax - dmin)/step)` bins over the data range; `Histogram[data, {min, max,
step}]` fixes the range and gives `ceil((max - min)/step)` bins. Counting: with
`bin_w = (bin_max - bin_min)/n_bins`, each value lands in bin
`(int)((v - bin_min)/bin_w)` clamped to `[0, n_bins-1]`, accumulated into an
`int` count array. Empty bins are skipped (no zero-height rectangle). Each
populated bin is drawn by the shared `emit_bar` as a 5-primitive bar over
`[x0, x0 + bin_w] x [0, count]`, coloured by `bar_color` indexed by the bin.
The result is an inert `Graphics[prims, opts...]` with `PlotRange` x = the bin
range and y = `[0, max_count*1.06 + 0.5]`.

**Data structures.** `ChartOpts`; a `calloc`'d `int counts[n_bins]`; `Expr**
prims` sized by the populated-bin count.

**Complexity / limits.** `O(n_vals + n_bins)`. Uniform-width bins only; Sturges'
count is hard-capped at 50; frequency counts only (no density/PDF
normalisation). Defaults `Axes -> True`, `AspectRatio -> 0.618` (`BarSpacing` is
parsed but unused in the histogram layout).
