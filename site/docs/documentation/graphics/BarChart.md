# BarChart

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`BarChart[{v1, v2, ..., vn}, opts...]`**

Draws a vertical bar chart: n bars at x = 1..n with heights v1..vn.

**`BarChart[{{v1,...}, {w1,...}, ...}, opts...]`**

Multiple grouped datasets, each in a distinct palette colour. Options: ChartStyle    color/style list cycling through bars (default: palette) ChartLabels   list of x-axis tick labels BarSpacing    gap fraction of bar width (default 0.2) Standard Graphics options (Axes, AspectRatio, Frame, PlotRange, PlotLabel, Background, ImageSize, …) pass through.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= BarChart[{3, 1, 4, 1, 5, 9, 2, 6}]
Out[1]= -Graphics-

In[2]:= BarChart[{3, -1, 4, -1, 5}]
Out[2]= -Graphics-
```

### Options (2)

```mathematica
In[3]:= BarChart[{2.5, 4.1, 3.3, 5.7}, ChartStyle -> {Red, Blue, Green, Orange}, ChartLabels -> {"Q1", "Q2", "Q3", "Q4"}]
Out[3]= -Graphics-

In[4]:= BarChart[{{1, 3, 2}, {4, 2, 5}}, BarSpacing -> 0.3]
Out[4]= -Graphics-
```

### Applications (5)

```mathematica
In[5]:= BarChart[{3, 1, 4, 1, 5, 9, 2, 6}]
Out[5]= -Graphics-

In[6]:= BarChart[{{1, 3, 2}, {4, 2, 5}}, BarSpacing -> 0.3]
Out[6]= -Graphics-

In[7]:= BarChart[{3, -1, 4, -1, 5}]
Out[7]= -Graphics-

In[8]:= BarChart[{2.5, 4.1, 3.3, 5.7}, ChartStyle -> {Red, Blue, Green, Orange}, ChartLabels -> {"Q1", "Q2", "Q3", "Q4"}]
Out[8]= -Graphics-

In[9]:= Length[Cases[BarChart[{3, 1, 4}], _Rectangle, Infinity]]
Out[9]= 3
```

## Algorithm

barchart.c — BarChart[data, opts...] and Histogram[data, opts...]

BarChart renders a vertical bar chart from explicit heights. Histogram bins numeric data and renders a frequency histogram. Both return Graphics[...] objects auto-displayed by the REPL.

BarChart[{v1,...,vn}, opts...]

```text
  n bars at x = 1..n with heights v1..vn, coloured via ChartStyle or
  the default palette.
```

BarChart[{{v1,...}, {w1,...}, ...}, opts...]

```text
  Multiple grouped datasets, each dataset in a distinct palette colour.
```

Histogram[data, opts...]

```text
Histogram[data, k, opts...]          k equal-width bins
Histogram[data, {step}, opts...]     bins of given width
```

Histogram[data, {min,max,step}, opts...] explicit range + width

Options (both):

```text
  ChartStyle   — color/style list cycling through bars
  ChartLabels  — label list for x-axis ticks
  BarSpacing   — gap fraction of bar width (default 0.2)
  PlotLabel, standard Graphics options pass through 
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [ChartStyle](../../other-advanced/ChartStyle/), [ChartLabels](../../other-advanced/ChartLabels/), [BarSpacing](../../other-advanced/BarSpacing/)

- Source: [`src/graphics/barchart.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/barchart.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)

## Notes & additional examples

### Notes

`BarChart` is `Protected` but **not** `HoldAll`, so its data is evaluated and each
element coerced to a real (non-numeric entries are dropped). A flat list gives `n`
bars at `x = 1..n`; a list of lists gives grouped sub-bars, one `palette_color` per
dataset. Each bar spans `[min(v, 0), max(v, 0)]`, so negatives draw downward.

`BarSpacing` (default 0.2) is the gap as a fraction of the group width; `ChartStyle`
cycles explicit colours; `ChartLabels` attaches category labels below each bar. Each
bar is drawn as a filled `Rectangle` plus a thin outline. Defaults are
`Axes -> True`, `AspectRatio -> 0.618`.
