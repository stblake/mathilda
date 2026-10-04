---
source: src/graphics/barchart.c
---
**Algorithm.** `builtin_barchart` is a plain `Protected` builtin (**not**
`HoldAll`): `args[0]` is evaluated and `extract_reals` coerces each element,
dropping non-numeric entries. A flat list is a single dataset of `n` bars centred
at `x = 1..n`; a list of lists is a grouped chart, one `palette_color` per
dataset. Layout: with `group_w = 1 - BarSpacing`, each group is `n_datasets`
sub-bars of width `group_w/n_datasets`, and a bar spans `[min(v, 0), max(v, 0)]`
so negative values draw downward from the axis. `bar_color` cycles an explicit
`ChartStyle` list (`i % len`), applies a single directive to all bars, or falls
back to the palette. `emit_bar` emits **five primitives per bar**: a fill colour
directive, a `Rectangle[{x0, y0}, {x1, y1}]`, a `GrayLevel[0.15]`, a
`Thickness[0.004]`, and a closed 5-point `Line` outline. `ChartLabels` adds a
`$BarChartLabels[{x, "label"}, ...]` metadata node for screen-space labels. The
result is an inert `Graphics[prims, opts...]` with `PlotRange` padded 5% in y.

**Data structures.** `ChartOpts` (from `split_chart_options`); per-dataset real
arrays from `extract_reals`; `Expr** prims` sized `n_datasets*n_bars*6 + ...`.

**Complexity / limits.** `O(total bars)`. Vertical bars only; a malformed inner
list in the multi-dataset form is replaced by a single zero bar. Defaults
`Axes -> True`, `AspectRatio -> 0.618`, `BarSpacing -> 0.2`.
