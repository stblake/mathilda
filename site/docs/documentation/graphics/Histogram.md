# Histogram

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Histogram[data, opts...]`**

Bins the numeric values in data and draws a frequency histogram. Bin count defaults to Sturges' rule: ceil(Log2\[n\]) + 1.

**`Histogram[data, k, opts...]`**

k equal-width bins.

**`Histogram[data, {step}, opts...]`**

Bins of width step.

**`Histogram[data, {min, max, step}, opts...]`**

Explicit range and width. Options: ChartStyle   color/style list cycling through bins BarSpacing   gap fraction of bin width (default 0.2) Standard Graphics options pass through.

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Histogram[Table[RandomReal[], {200}]]
Out[1]= -Graphics-

In[2]:= data = Table[RandomReal[], {200}]
Out[2]= {0.898168, 0.822246, 0.486932, 0.908284, 0.428175, 0.932029, 0.131627, 0.853968, 0.0497704, 0.584082, 0.302085, 0.137246, 0.701554, 0.55101, 0.282818, 0.717397, 0.420001, 0.248719, 0.352146, 0.834998, 0.0979758, 0.104027, 0.709625, 0.835462, 0.729187, 0.597183, 0.53771, 0.946709, 0.41, 0.305271, 0.0779037, 0.334263, 0.354373, 0.00963086, 0.591771, 0.523084, 0.0360014, 0.392808, 0.183878, 0.223776, 0.23013, 0.681739, 0.636806, 0.194441, 0.203807, 0.725877, 0.375132, 0.867256, 0.431451, 0.810405, 0.772533, 0.4733, 0.608341, 0.548074, 0.650504, 0.767579, 0.402964, 0.903085, 0.255249, 0.488619, 0.215692, 0.483165, 0.0850616, 0.975261, 0.664815, 0.698508, 0.242494, 0.177173, 0.152351, 0.678393, 0.764858, 0.0119961, 0.724651, 0.053345, 0.38038, 0.546757, 0.101937, 0.691765, 0.425465, 0.10441, 0.236588, 0.629644, 0.256469, 0.359752, 0.243209, 0.765368, 0.0625062, 0.395556, 0.463383, 0.367876, 0.242292, 0.280831, 0.222355, 0.95352, 0.916981, 0.291277, 0.790124, 0.656521, 0.496403, 0.36685, 0.605281, 0.190074, 0.368757, 0.750386, 0.218605, 0.854503, 0.0631354, 0.492749, 0.78421, 0.067381, 0.404588, 0.887172, 0.986161, 0.109416, 0.0188128, 0.857345, 0.199804, 0.397022, 0.964901, 0.137297, 0.545313, 0.550262, 0.62597, 0.508226, 0.177263, 0.900855, 0.670446, 0.924974, 0.986038, 0.457149, 0.730116, 0.357833, 0.264861, 0.489474, 0.0695345, 0.228527, 0.0651126, 0.878594, 0.0886117, 0.858396, 0.870913, 0.939162, 0.260281, 0.0733191, 0.756655, 0.905309, 0.0368469, 0.89336, 0.385519, 0.597237, 0.888194, 0.429421, 0.41403, 0.525518, 0.890559, 0.657609, 0.197113, 0.833592, 0.144836, 0.16589, 0.738772, 0.095379, 0.6518, 0.263895, 0.144144, 0.276424, 0.933648, 0.342755, 0.365329, 0.640436, 0.346002, 0.595686, 0.0739085, 0.883215, 0.464774, 0.628304, 0.849411, 0.490909, 0.376777, 0.536012, 0.539017, 0.950868, 0.83077, 0.207706, 0.553058, 0.753565, 0.477569, 0.333518, 0.269935, 0.448931, 0.887562, 0.130198, 0.396262, 0.902854, 0.460089, 0.457556, 0.342479, 0.952447, 0.504037, 0.291491}

In[3]:= Histogram[data, 20]
Out[3]= -Graphics-

In[4]:= Histogram[data, {0.1}]
Out[4]= -Graphics-

In[5]:= Histogram[data, {0, 1, 0.05}]
Out[5]= -Graphics-
```

### Applications (5)

```mathematica
In[6]:= Histogram[{1, 2, 2, 3, 3, 3, 4, 4, 4, 4}]
Out[6]= -Graphics-

In[7]:= Histogram[Range[100], 10]
Out[7]= -Graphics-

In[8]:= Histogram[Range[100], {25}]
Out[8]= -Graphics-

In[9]:= Histogram[Range[100], {0, 100, 25}]
Out[9]= -Graphics-

In[10]:= Length[Cases[Histogram[Range[100], 10], _Rectangle, Infinity]]
Out[10]= 10
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

**Attributes:** `Protected`.

## References

**See also:** [ChartStyle](../../other-advanced/ChartStyle/), [BarSpacing](../../other-advanced/BarSpacing/)

- H. A. Sturges, *The Choice of a Class Interval*, J. Amer. Statist. Assoc. **21** (1926) 65-66 — the default bin count.
- Source: [`src/graphics/barchart.c`](https://github.com/stblake/mathilda/blob/main/src/graphics/barchart.c)
- Specification: [`docs/spec/builtins/graphics.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphics.md)

## Notes & additional examples

### Notes

`Histogram` lives alongside `BarChart` and shares its bar-drawing machinery; it is
`Protected` but **not** `HoldAll`. The bin specification is flexible:
`Histogram[data]` uses **Sturges' rule** (`ceil(log2 n) + 1` bins, clamped to
`[2, 50]`); `Histogram[data, k]` gives `k` bins; `Histogram[data, {step}]` fixes the
bin width; `Histogram[data, {min, max, step}]` fixes the range and width.

Each value is dropped into its bin and the bin counts become bar heights; empty
bins are skipped (no zero-height rectangle). Only equal-width frequency bins are
supported — there is no density/PDF normalisation. Defaults are `Axes -> True`,
`AspectRatio -> 0.618`.
